#include "flipsolver2d.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <limits>
#include <memory>
#include <ostream>
#include <queue>

#include "PressureIPPCoeficients.h"
#include "grid2d.h"
#include "index2d.h"
#include "linearindexable2d.h"
#include "linearsolver.h"
#include "markerparticlesystem.h"
#include "materialgrid.h"

#include "mathfuncs.h"
#include "pressuredata.h"
#include "threadpool.h"
#include "viscositymodel.h"

#include <Eigen/IterativeLinearSolvers>

FlipSolver::FlipSolver(const FlipSolverParameters *p) :
    LinearIndexable2d(p->gridSizeI, p->gridSizeJ),
    m_frameNumber(0),
    m_markerParticles(p->gridSizeI, p->gridSizeJ,3),
    m_fluidVelocityGrid(p->gridSizeI, p->gridSizeJ),
    m_savedFluidVelocityGrid(p->gridSizeI, p->gridSizeJ),
    m_materialGrid(p->gridSizeI,p->gridSizeJ, FluidMaterial::SINK),
    m_solidSdf(p->gridSizeI,p->gridSizeJ),
    m_fluidSdf(p->gridSizeI, p->gridSizeJ),
    m_knownCenteredParams(p->gridSizeI,p->gridSizeJ, false, OOBStrategy::OOB_CONST, true),
    m_viscosityGrid(p->gridSizeI,p->gridSizeJ,0.f,OOBStrategy::OOB_EXTEND, 0.f, Vec3(0.5f,0.5f)),
    m_emitterId(p->gridSizeI, p->gridSizeJ, -1),
    m_solidId(p->gridSizeI, p->gridSizeJ,-1),
    m_fluidParticleCounts(p->gridSizeI, p->gridSizeJ),
    m_divergenceControl(p->gridSizeI,p->gridSizeJ, 0.f, OOBStrategy::OOB_CONST, 0.f),
    m_densityGrid(p->gridSizeI,p->gridSizeJ, 0.f),
    m_testGrid(p->gridSizeI,p->gridSizeJ),
    m_rhs(),
    m_solverResult(),
    m_stepDt(1.f / p->fps),
    m_frameDt(1.f / p->fps), m_dx(p->dx),
    m_fluidDensity(p->fluidDensity),
    m_seed(p->seed),
    m_particlesPerCell(p->particlesPerCell),
    m_globalAcceleration(p->globalAcceleration),
    m_resolution(p->resolution),
    m_fps(p->fps),
    m_maxSubsteps(p->maxSubsteps),
    m_picRatio(p->picRatio),
    m_cflNumber(p->cflNumber),
    m_particleScale(p->particleScale),
    m_pcgIterLimit(p->pcgIterLimit),
    m_domainSizeI(p->domainSizeI),
    m_domainSizeJ(p->domainSizeJ),
    m_sceneScale(p->sceneScale),
    m_projectTolerance(1e-6),
    m_viscosityEnabled(p->viscosityEnabled),
    m_simulationMethod(p->simulationMethod),
    m_parameterHandlingMethod(p->parameterHandlingMethod),
    m_pressureMatrix(0,*this,0.0),
    m_pressurePrecond(0,*this)
{
    Eigen::initParallel();
    Eigen::setNbThreads(ThreadPool::i()->threadCount());
    m_randEngine = std::mt19937(p->seed);
    m_testValuePropertyIndex = m_markerParticles.addParticleProperty<float>();
    m_projectTolerance = m_viscosityEnabled? 1e-6 : 1e-2;
    m_rhs.resize(linearSize());
    m_solverResult.resize(linearSize());

    if(p->useHeavyViscosity)
    {
        m_viscosityModel = std::make_shared<HeavyViscosityModel>();
        std::cout << "using heavy model" << std::endl;
    }
    else
    {
        m_viscosityModel = std::make_shared<LightViscosityModel>();
    }
}

FlipSolver::~FlipSolver()
{
    ThreadPool::i()->wait();
}

void FlipSolver::project()
{
    //using precond = Eigen::IncompleteLUT<double>;
    std::vector<double> rhs(cellCount(), 0.0);
    std::vector<double> solverResult(cellCount(), 0.0);

    calcPressureRhs(rhs);
//    //debug() << "Calculated rhs: " << rhs;
    //Eigen::BiCGSTAB<Eigen::SparseMatrix<double>,precond> solver;

    if(!m_pressureSolver.solve(m_pressureMatrix,m_pressurePrecond,
                                solverResult,rhs,
                                m_pcgIterLimit, m_projectTolerance)) {
        std::cout << "Pressure solver solving failed! Expect bogus pressures\n";

    }

//    if(!m_pcgSolver.mfcgSolve(provider,m_pressures.data(),m_rhs,m_pcgIterLimit,1e-2))
//    {
//        std::cout << "PCG Solver pressure: Iteration limit exhaustion!\n";
//    }

//    if(anyNanInf(m_pressures.data()))
//    {
//        std::cout << "NaN or inf in pressures vector!\n" << std::flush;
//    }

    //debug() << "pressures = " << pressures;

    applyPressuresToVelocityField(solverResult);
}

void FlipSolver::applyViscosity()
{
    m_viscosityModel->apply(m_fluidVelocityGrid,
                            m_viscosityGrid,
                            m_materialGrid,
                            m_stepDt,
                            m_dx,
                            m_fluidDensity);
}

void FlipSolver::advect()
{
    std::vector<Range> ranges = ThreadPool::i()->splitRange(m_markerParticles.bins().linearSize());
    std::vector<std::vector<RebinRecord>> rebinSets(m_markerParticles.bins().linearSize());

    for(size_t i = 0; i < ranges.size(); i++)
    {
        ThreadPool::i()->enqueue(&FlipSolver::advectThread,this,ranges.at(i),std::ref(rebinSets));
    }
    ThreadPool::i()->wait();

    for(ssize_t binIdx = 0; binIdx < m_markerParticles.bins().linearSize(); binIdx++)
    {
        for(RebinRecord r : rebinSets[binIdx])
        {
            m_markerParticles.scheduleRebin(binIdx,r);
        }
    }

    if(m_parameterHandlingMethod == GRID)
    {
        eulerAdvectParameters();
    }
    //std::cout << "Advection done in max " << maxSubsteps << " substeps" << std::endl;
}

void FlipSolver::densityCorrection()
{
    updateDensityGrid();

    m_pressureRhs.resize(cellCount());
    m_pressureSolverResult.resize(cellCount());

    calcDensityCorrectionRhs(m_pressureRhs);

    if(!m_pressureSolver.solve(m_pressureMatrix, m_pressurePrecond,
                               m_pressureSolverResult, m_pressureRhs, m_pcgIterLimit,1e-4))
    {
        std::cout << "Density solver solving failed!\n";
        return;
    }

    adjustParticlesByDensity();
    std::cout << "Density correction done\n";
}

void FlipSolver::updateDensityGrid()
{
    Grid2d<float> centeredWeights(m_sizeI,m_sizeJ,1e-10f);
    m_densityGrid.fill(0.f);

    std::vector<Range> ranges = ThreadPool::i()->splitRange(m_densityGrid.linearSize(),1,8);
    for(Range& range : ranges)
    {
        ThreadPool::i()->enqueue(&FlipSolver::updateDensityGridThread,this,range,std::ref(centeredWeights));
    }
    ThreadPool::i()->wait();
}

void FlipSolver::updateDensityGridThread(Range r, Grid2d<float> &centeredWeights)
{
    const float cellVolume = m_dx * m_dx * m_dx;
    const float particleMass = (m_fluidDensity * cellVolume) / (static_cast<float>(m_particlesPerCell));
    auto& bins = m_markerParticles.bins();

    for(ssize_t idx = r.start; idx < r.end; idx++)
    {
        Index2d cellIdx = m_densityGrid.index2d(idx);
        for(ssize_t binIdx : m_markerParticles.binsForGridCell(cellIdx))
        {
            if(binIdx == -1)
            {
                continue;
            }

            ParticleBin& currentBin = m_markerParticles.bins().data()[binIdx];

            for(size_t particleIdx = 0; particleIdx < currentBin.size(); particleIdx++)
            {
                if(currentBin.markedForDeath(particleIdx))
                {
                    continue;
                }

                Vec3 p = currentBin.particlePosition(particleIdx);
                float weightCentered = simmath::bilinearHat(p.x() - static_cast<float>(cellIdx.i) - 0.5f,
                                                            p.y() - static_cast<float>(cellIdx.j) - 0.5f);
                if(std::abs(weightCentered) > 1e-6f)
                {
                    centeredWeights.at(cellIdx) += weightCentered;
                    m_densityGrid.at(cellIdx) += weightCentered * particleMass;
                }
            }
        }
        m_densityGrid.at(cellIdx) /= cellVolume;
        if(m_materialGrid.isFluid(cellIdx) &&
            (m_materialGrid.isEmpty(cellIdx.i+1,cellIdx.j) ||
             m_materialGrid.isEmpty(cellIdx.i-1,cellIdx.j) ||
             m_materialGrid.isEmpty(cellIdx.i,cellIdx.j+1) ||
             m_materialGrid.isEmpty(cellIdx.i,cellIdx.j-1)))
        {
            m_densityGrid.at(cellIdx) = std::clamp(m_densityGrid.at(cellIdx),
                                                static_cast<float>(m_fluidDensity),
                                                std::numeric_limits<float>::max());
        }
        //m_testGrid.at(cellIdx) = m_densityGrid.at(cellIdx);
    }
}

void FlipSolver::adjustParticlesByDensity()
{
    std::vector<Range> ranges = ThreadPool::i()->splitRange(m_markerParticles.bins().linearSize());
    for(Range& range : ranges)
    {
        ThreadPool::i()->enqueue(&FlipSolver::adjustParticlesByDensityThread,this,range);
    }
    ThreadPool::i()->wait();
}

void FlipSolver::adjustParticlesByDensityThread(Range r)
{
    const float scale = (m_stepDt * m_stepDt) / (m_fluidDensity * m_dx * m_dx);
    for(size_t binIdx = r.start; binIdx < r.end; binIdx++)
    {
        ParticleBin& currentBin = m_markerParticles.bins().data().at(binIdx);

        for(size_t particleIdx = 0; particleIdx < currentBin.size(); particleIdx++)
        {
            if(currentBin.markedForDeath(particleIdx))
            {
                continue;
            }

            Vec3& position = currentBin.particlePosition(particleIdx);
            ssize_t iCorr = std::clamp<ssize_t>(static_cast<int>(position.x() - 0.5f),0,m_sizeI);
            ssize_t jCorr = std::clamp<ssize_t>(static_cast<int>(position.y() - 0.5f),0,m_sizeJ);
            ssize_t i = std::clamp<ssize_t>(static_cast<int>(position.x()),0,m_sizeI);
            ssize_t j = std::clamp<ssize_t>(static_cast<int>(position.y()),0,m_sizeJ);

            float pCurrentI = m_pressureSolverResult[linearIndex(iCorr,j)];
            float pCurrentJ = m_pressureSolverResult[linearIndex(i,jCorr)];

            float pI = m_pressureSolverResult[linearIndex(iCorr+1,j)];
            float pJ = m_pressureSolverResult[linearIndex(i,jCorr+1)];
            if(m_materialGrid.isSolid(iCorr+1,j) || m_materialGrid.isSolid(iCorr,j))
            {
                pI = pCurrentI;
            }

            if(m_materialGrid.isSolid(i,jCorr+1) || m_materialGrid.isSolid(i,jCorr))
            {
                pJ = pCurrentJ;
            }

            float pgradU = pI - pCurrentI;
            float pgradV = pJ - pCurrentJ;

            position.x() += pgradU*scale;
            position.y() += pgradV*scale;
        }
    }
}

void FlipSolver::advectThread(Range range, std::vector<std::vector<RebinRecord>> &rebinningSets)
{
    //std::vector<float>& testValues = m_markerParticles.particleProperties<float>(m_testValuePropertyIndex);

    for(ssize_t binIdx = range.start; binIdx < range.end; binIdx++)
    {
        std::vector<RebinRecord>& rebinningSet = rebinningSets.at(binIdx);
        std::vector<Vec3>& positions = m_markerParticles.binForBinIdx(binIdx).positions();

        for(ssize_t particleIdx = 0; particleIdx < positions.size(); particleIdx++)
        {
            Vec3& position = positions[particleIdx];
            const ssize_t oldBinIdx = m_markerParticles.gridToBinIdx(position);
            position = rk4Integrate(position, m_fluidVelocityGrid, m_stepDt);
            if(m_solidSdf.interpolateAt(position.x(),position.y()) < 0.f)
            {
                position = m_solidSdf.closestSurfacePoint(position);
            }
            ssize_t pI = simmath::integr(position.x());
            ssize_t pJ = simmath::integr(position.y());
            //testValues[i] = false;
            const ssize_t newBinIdx = m_markerParticles.gridToBinIdx(position);
            if(!inBounds(pI,pJ) || m_materialGrid.isSink(pI,pJ))
            {
                m_markerParticles.markForDeath(binIdx,particleIdx);
            }
            else if(newBinIdx != oldBinIdx)
            {
                rebinningSet.push_back(RebinRecord(particleIdx,newBinIdx));
                //testValues[i] = true;
            }
        }
    }
}

void FlipSolver::eulerAdvectionThread(Range range, const Grid2d<float> &inputGrid, Grid2d<float> &outputGrid)
{
    std::vector<float>& dataOut = outputGrid.data();
    std::vector<int>& emitterIdx = m_emitterId.data();
    for(ssize_t idx = range.start; idx < range.end; idx++)
    {
        Index2d i2d = outputGrid.index2d(idx);
        Vec3 currentPos = Vec3(i2d.i, i2d.j);
        Vec3 prevPos = rk4Integrate(currentPos,m_fluidVelocityGrid, -m_stepDt);
        dataOut[idx] = inputGrid.interpolateAt(prevPos);
    }
}

void FlipSolver::pruneParticles()
{
    for(ParticleBin& bin : m_markerParticles.bins().data())
    {
        bin.pruneParticles();
    }
}

void FlipSolver::particleUpdate()
{
    Grid2d<float>& prevU = m_savedFluidVelocityGrid.velocityGridU();
    Grid2d<float>& prevV = m_savedFluidVelocityGrid.velocityGridV();
    for(size_t binIdx = 0; binIdx < m_markerParticles.bins().linearSize(); binIdx++)
    {
        ParticleBin& currentBin = m_markerParticles.bins().data().at(binIdx);

        for(size_t particleIdx = 0; particleIdx < currentBin.size(); particleIdx++)
        {
            Vec3 &position = currentBin.particlePosition(particleIdx);
            Vec3 &velocity = currentBin.particleVelocity(particleIdx);

            Vec3 oldVelocity(prevU.interpolateAt(position.x(),position.y()),
                               prevV.interpolateAt(position.x(),position.y()));
            Vec3 newVelocity = m_fluidVelocityGrid.velocityAt(position) ;
    //        if(oldVelocity.distFromZero() > (SimSettings::cflNumber() / m_stepDt))
    //        {
    //            velocity = newVelocity;
    //        }
    //        else
    //        {
                velocity = m_picRatio * newVelocity +
                        (1.f-m_picRatio) * (velocity + newVelocity - oldVelocity);
    //        }
        }
    }
}

void FlipSolver::afterTransfer()
{
    for (ssize_t i = 0; i < m_sizeI; i++)
    {
        for (ssize_t j = 0; j < m_sizeJ; j++)
        {
            if(m_materialGrid.isSource(i,j))
            {
                int emitterId = m_emitterId.at(i,j);
                m_viscosityGrid.at(i,j) = m_sources[emitterId].viscosity();
               if(m_sources[emitterId].velocityTransfer())
               {
                    m_fluidVelocityGrid.setU(i,j,m_sources[emitterId].velocity().x() / m_dx);
                    m_fluidVelocityGrid.setV(i,j,m_sources[emitterId].velocity().y() / m_dx);
                    m_fluidVelocityGrid.setUValidity(i,j,true);
                    m_fluidVelocityGrid.setVValidity(i,j,true);
               }
            }
        }
    }
}

void FlipSolver::step()
{
    advect();
    m_stats.endStage(ADVECTION);

    m_pressureMatrix = getPressureProjectionMatrix();
    m_pressurePrecond = getIPPCoefficients(m_pressureMatrix);
    m_stats.endStage(DECOMPOSITION);

    pruneParticles();
    m_markerParticles.rebinParticles();
    m_stats.endStage(PARTICLE_REBIN);

    if(!m_viscosityEnabled)
    {
        //Assumption: particles are adjusted not enough to require rebinning
        densityCorrection();
        m_stats.endStage(DENSITY);
    }

    gridUpdate();
    m_stats.endStage(GRID_UPDATE);

    afterTransfer();
    extrapolateLevelsetInside(m_fluidSdf);
    m_fluidVelocityGrid.extrapolate(10);

    m_savedFluidVelocityGrid = m_fluidVelocityGrid;
    applyBodyForces();
    m_stats.endStage(AFTER_TRANSFER);

    //m_pressureSolver.setTolerance(1e-6);
    project();
    m_stats.endStage(PRESSURE);

    updateVelocityFromSolids();
    if(m_viscosityEnabled)
    {
        applyViscosity();
        m_stats.endStage(VISCOSITY);
        project();
        m_stats.endStage(REPRESSURE);
    }
    m_fluidVelocityGrid.extrapolate(10);
    particleUpdate();
    m_stats.endStage(PARTICLE_UPDATE);

    countParticles();
    reseedParticles();
    m_stats.endStage(PARTICLE_RESEED);
}

void FlipSolver::stepFrame()
{
    if(m_frameNumber == 0)
    {
        firstFrameInit();
    }
    m_testGrid.fill(0.f);
    float substepTime = 0.f;
    bool finished = false;
    int substepCount = 0;
    m_stats.reset();

    while(!finished)
    {
        float vel = maxParticleVelocity();
        float maxSubstepSize = m_cflNumber/(vel + 1e-15f);
        if(substepTime + maxSubstepSize >= m_frameDt ||
                substepCount == (m_maxSubsteps - 1))
        {
            maxSubstepSize = m_frameDt - substepTime;
            finished = true;
        }
        else if(substepTime + 2.f*maxSubstepSize >= m_frameDt)
        {
            maxSubstepSize = 0.5f*(m_frameDt - substepTime);
        }
        m_stepDt = maxSubstepSize;
        std::cout << "Substep " << substepCount << " substep dt: " << m_stepDt << " vel " << vel << std::endl;
        step();
        m_stats.addSubstep();
        substepTime += maxSubstepSize;
        substepCount++;
        if(substepCount > 50) break;
    }
    m_stats.endFrame();
    m_frameNumber++;
}

void FlipSolver::updateSolids()
{
    for (size_t i = 0; i < m_sizeI; i++)
    {
        for (size_t j = 0; j < m_sizeJ; j++)
        {
            float dx = m_dx;
            float dist = std::numeric_limits<float>::max();
            int minIdx = 0;
            for(int solidIdx = 0; solidIdx < m_obstacles.size(); solidIdx++)
            {
                float sdf = m_obstacles[solidIdx].geometry().signedDistance(
                            (static_cast<float>(i)+0.5)*dx,
                            (static_cast<float>(j)+0.5)*dx) / dx;
                if(sdf < dist)
                {
                    minIdx = solidIdx;
                    dist = sdf;
                }
            }
            m_solidSdf.at(i,j) = dist;
            if(dist < 0)
            {
                m_materialGrid.at(i,j) = FluidMaterial::SOLID;
                m_solidId.at(i,j) = minIdx;
            }
        }
    }
}

void FlipSolver::updateSources()
{
    const float dx = m_dx;
    const float hdx = dx / 2.0;
    for (size_t i = 0; i < m_sizeI; i++)
    {
        for (size_t j = 0; j < m_sizeJ; j++)
        {
            for(int emitterIdx = 0; emitterIdx < m_sources.size(); emitterIdx++)
            {
                Emitter& e = m_sources[emitterIdx];
                if(e.geometry().signedDistance(i*dx + hdx,j*dx + hdx) <= 0.f)
                {
                    m_materialGrid.at(i,j) = FluidMaterial::SOURCE;
                    m_emitterId.at(i,j) = emitterIdx;
                    m_divergenceControl.at(i,j) = e.divergence();
                }
            }
        }
    }
}

void FlipSolver::updateSinks()
{
    float dx = m_dx;
    for (size_t i = 0; i < m_sizeI; i++)
    {
        for (size_t j = 0; j < m_sizeJ; j++)
        {
            for(Sink& s : m_sinks)
            {
                if(s.geo().signedDistance((static_cast<float>(i)+0.5f)*dx,(static_cast<float>(j)+0.5f)*dx) <= 0.f)
                {
                    m_materialGrid.at(i,j) = FluidMaterial::SINK;
                    m_divergenceControl.at(i,j) = s.divergence();
                }
            }
        }
    }
}

void FlipSolver::updateInitialFluid()
{
    float dx = m_dx;
    for (size_t i = 0; i < m_sizeI; i++)
    {
        for (size_t j = 0; j < m_sizeJ; j++)
        {
            for(Emitter& e : m_initialFluid)
            {
                if(e.geometry().signedDistance((static_cast<float>(i)+0.5f)*dx,(static_cast<float>(j)+0.5f)*dx) <= 0.f)
                {
                    m_materialGrid.at(i,j) = FluidMaterial::FLUID;
                    m_viscosityGrid.at(i,j) = e.viscosity();
                }
            }
        }
    }
}

size_t FlipSolver::gridSizeI()
{
    return m_sizeI;
}

size_t FlipSolver::gridSizeJ()
{
    return m_sizeJ;
}

void FlipSolver::addGeometry(Obstacle& geometry)
{
    m_obstacles.push_back(geometry);
}

void FlipSolver::addSource(Emitter &emitter)
{
    m_sources.push_back(emitter);
}

void FlipSolver::addSink(Sink &geometry)
{
    m_sinks.push_back(geometry);
}

void FlipSolver::addInitialFluid(Emitter &emitter)
{
    m_initialFluid.push_back(emitter);
}

int FlipSolver::frameNumber()
{
    return m_frameNumber;
}

void FlipSolver::reseedParticles()
{
    // std::vector<float>& particleViscosities = m_markerParticles.particleProperties<float>
    //                                           (m_viscosityPropertyIndex);

    for (size_t i = 0; i < m_sizeI; i++)
    {
        for (size_t j = 0; j < m_sizeJ; j++)
        {
            int particleCount = m_fluidParticleCounts.at(i,j);
            if(particleCount > m_particlesPerCell*2)
            {
                std::cout << "too many particles " << particleCount << " at " << i << ' ' << j;
            }
            //std::cout << particleCount << " at " << i << " , " << j << std::endl;
            int additionalParticles = (m_particlesPerCell / 2) - particleCount;
            if(additionalParticles <= 0)
            {
                continue;
            }

            ParticleBin& bin = m_markerParticles.binForGridIdx(i,j);
            auto& viscVec = bin.particleProperties<float>(m_viscosityPropertyIndex);

            if(m_materialGrid.isSource(i,j))
            {
                int emitterId = m_emitterId.at(i,j);
                for(int p = 0; p < additionalParticles; p++)
                {
                    Vec3 pos = jitteredPosInCell(i,j);
                    Vec3 velocity = Vec3();
                    if(m_sources[emitterId].velocityTransfer())
                    {
                        velocity = m_fluidVelocityGrid.velocityAt(pos);
                    }
                    float viscosity = m_sources[emitterId].viscosity();
                    size_t pIdx = m_markerParticles.binForGridPosition(pos).addMarkerParticle(pos,velocity);
                    //size_t pIdx = m_markerParticles.addMarkerParticle(pos,velocity);
                    viscVec[pIdx] = viscosity;
                }
            }
//            else if(m_fluidSdf.at(i,j) < -1.f)
//            {
//                for(int p = 0; p < additionalParticles; p++)
//                {
//                    Vertex pos = jitteredPosInCell(i,j);
//                    Vertex velocity = m_fluidVelocityGrid.velocityAt(pos);
//                    float viscosity = m_viscosityGrid.interpolateAt(pos);
//                    addMarkerParticle(MarkerParticle{pos,velocity,viscosity});
//                }
//            }
        }
    }
}

void FlipSolver::seedInitialFluid()
{
    // std::vector<float>& particleViscosities = m_markerParticles.particleProperties<float>
    //                                           (m_viscosityPropertyIndex);
    for (size_t i = 0; i < m_sizeI; i++)
    {
        for (size_t j = 0; j < m_sizeJ; j++)
        {
            if(m_materialGrid.isStrictFluid(i,j))
            {
                ParticleBin& bin = m_markerParticles.binForGridIdx(i,j);
                auto& viscVec = bin.particleProperties<float>(m_viscosityPropertyIndex);

                for(int p = 0; p < m_particlesPerCell; p++)
                {
                    Vec3 pos = jitteredPosInCell(i,j);
                    Vec3 velocity = m_fluidVelocityGrid.velocityAt(pos);
                    float viscosity = m_viscosityGrid.interpolateAt(pos);

                    size_t pIdx = m_markerParticles.binForGridPosition(pos).addMarkerParticle(pos,velocity);
                    viscVec.at(pIdx) = viscosity;
                }
            }
        }
    }
}

std::vector<Obstacle> &FlipSolver::geometryObjects()
{
    return m_obstacles;
}

std::vector<Emitter> &FlipSolver::sourceObjects()
{
    return m_sources;
}

std::vector<Sink> &FlipSolver::sinkObjects()
{
    return m_sinks;
}

MarkerParticleSystem &FlipSolver::markerParticles()
{
    return m_markerParticles;
}

const MaterialGrid &FlipSolver::materialGrid() const
{
    return m_materialGrid;
}

const StaggeredVelocityGrid &FlipSolver::fluidVelocityGrid() const
{
    return m_fluidVelocityGrid;
}

const SdfGrid &FlipSolver::fluidSdf() const
{
    return m_fluidSdf;
}

const SdfGrid &FlipSolver::solidSdf() const
{
    return m_solidSdf;
}

const Grid2d<float> &FlipSolver::testGrid() const
{
    return m_testGrid;
}

const SolverTimeStats &FlipSolver::timeStats() const
{
    return m_stats;
}

double FlipSolver::divergenceAt(size_t i, size_t j)
{
    return static_cast<double>(m_fluidVelocityGrid.u(i+1,j) - m_fluidVelocityGrid.u(i,j)
                               + m_fluidVelocityGrid.v(i,j+1) - m_fluidVelocityGrid.v(i,j)
                               + m_divergenceControl.at(i,j));
}

std::vector<int> FlipSolver::validSolidNeighborIds(ssize_t i, ssize_t j)
{
    std::vector<int> output;
    output.reserve(8);
    for (ssize_t iOffset = -1; iOffset < 2; iOffset++)
    {
        for (ssize_t jOffset = -1; jOffset < 2; jOffset++)
        {
            if(m_solidId.inBounds(i + iOffset, j + jOffset))
            {
                int id = m_solidId.at(i + iOffset, j + jOffset);
                if(id != -1)
                {
                    output.push_back(m_solidId.at(i + iOffset, j + jOffset));
                }
            }
        }
    }

    return output;
}

void FlipSolver::firstFrameInit()
{
    updateInitialFluid();
    seedInitialFluid();
}

IndexedPressureParameters FlipSolver::getPressureProjectionMatrix()
{
    const double scale = m_stepDt / (m_fluidDensity * m_dx * m_dx);

    IndexedPressureParameters output(linearSize()*0.33, *this, scale);

    LinearIndexable2d& indexer = *dynamic_cast<LinearIndexable2d*>(this);

    std::vector<Range> threadRanges = ThreadPool::i()->splitRange(linearSize());
    size_t currRangeIdx = 0;

    for(size_t i = 0; i < m_sizeI; i++)
    {
        for(size_t j = 0; j < m_sizeJ; j++)
        {
            const size_t linIdx = indexer.linearIndex(i,j);

            if(!m_materialGrid.isFluid(i,j))
            {
                if(linIdx >= threadRanges.at(currRangeIdx).end)
                {
                    output.endThreadDataRange();
                    currRangeIdx++;
                }
                continue;
            }

            IndexedPressureParameterUnit unit;
            unit.unitIndex = linIdx;

            const ssize_t linIdxPosAx = indexer.linearIdxOfOffset(linIdx,1,0);
            const ssize_t linIdxPosAy = indexer.linearIdxOfOffset(linIdx,0,1);
            const ssize_t linIdxNegAx = indexer.linearIdxOfOffset(linIdx,-1,0);
            const ssize_t linIdxNegAy = indexer.linearIdxOfOffset(linIdx,0,-1);

            double diag = 0.0;
            //X Neighbors
            if(m_materialGrid.isFluid(i-1,j))
            {
                unit.nonsolidNeighborCount++;
                unit.setNeighbor(I_NEG_NEIGHBOR_BIT, inBounds(linIdxNegAx));
            }else if(m_materialGrid.isEmpty(i-1,j))
            {
                unit.nonsolidNeighborCount++;
            }

            if(m_materialGrid.isFluid(i+1,j))
            {
                unit.nonsolidNeighborCount++;
                unit.setNeighbor(I_POS_NEIGHBOR_BIT, inBounds(linIdxPosAx));
            } else if(m_materialGrid.isEmpty(i+1,j))
            {
                unit.nonsolidNeighborCount++;
            }

            //Y Neighbors
            if(m_materialGrid.isFluid(i,j-1))
            {
                unit.nonsolidNeighborCount++;
                unit.setNeighbor(J_NEG_NEIGHBOR_BIT, inBounds(linIdxNegAy));
            }else if(m_materialGrid.isEmpty(i,j-1))
            {
                unit.nonsolidNeighborCount++;
            }

            if(m_materialGrid.isFluid(i,j+1))
            {
                unit.nonsolidNeighborCount++;
                unit.setNeighbor(J_POS_NEIGHBOR_BIT, inBounds(linIdxPosAy));
            } else if(m_materialGrid.isEmpty(i,j+1))
            {
                unit.nonsolidNeighborCount++;
            }

            if(linIdx >= threadRanges.at(currRangeIdx).end)
            {
                output.endThreadDataRange();
                currRangeIdx++;
            }

            output.add(unit);
        }
    }

    if(currRangeIdx < threadRanges.size())
    {
        output.endThreadDataRange();
    }

    return output;
}

IndexedIPPCoefficients FlipSolver::getIPPCoefficients(const IndexedPressureParameters& mat)
{
    const double scale = m_stepDt / (m_fluidDensity * m_dx * m_dx);

    IndexedIPPCoefficients output(linearSize()*0.33, *this);

    LinearIndexable2d& indexer = *dynamic_cast<LinearIndexable2d*>(this);

    std::vector<Range> threadRanges = ThreadPool::i()->splitRange(linearSize());
    size_t currRangeIdx = 0;

    size_t matEntryIdx = 0;

    for(size_t i = 0; i < m_sizeI; i++)
    {
        for(size_t j = 0; j < m_sizeJ; j++)
        {
            const size_t linIdx = indexer.linearIndex(i,j);

            if(!m_materialGrid.isFluid(i,j))
            {
                if(linIdx >= threadRanges.at(currRangeIdx).end)
                {
                    output.endThreadDataRange();
                    currRangeIdx++;
                }
                continue;
            }

            IndexedIPPCoefficientUnit unit;
            unit.unitIndex = linIdx;

            const double invscale = scale/(mat.data().at(matEntryIdx).nonsolidNeighborCount*scale);

            const ssize_t iNegLinIdx = indexer.linearIdxOfOffset(linIdx,-1,0);
            const ssize_t iPosLinIdx = indexer.linearIdxOfOffset(linIdx,1,0);
            const ssize_t jNegLinIdx = indexer.linearIdxOfOffset(linIdx,0,-1);
            const ssize_t jPosLinIdx = indexer.linearIdxOfOffset(linIdx,0,1);

            unit.iNeg = 1.0/m_materialGrid.nonsolidNeighborCount(iNegLinIdx);
            unit.iPos = 1.0/m_materialGrid.nonsolidNeighborCount(iPosLinIdx);
            unit.jNeg = 1.0/m_materialGrid.nonsolidNeighborCount(jNegLinIdx);
            unit.jPos = 1.0/m_materialGrid.nonsolidNeighborCount(jPosLinIdx);

            if(linIdx >= threadRanges.at(currRangeIdx).end)
            {
                output.endThreadDataRange();
                currRangeIdx++;
            }

            output.add(unit);
            matEntryIdx++;
        }
    }

    if(currRangeIdx != threadRanges.size())
    {
        output.endThreadDataRange();
    }

    return output;
}

size_t FlipSolver::particleCount()
{
    return m_markerParticles.particleCount();
}

size_t FlipSolver::cellCount()
{
    return m_sizeI * m_sizeJ;
}

void FlipSolver::calcPressureRhs(std::vector<double> &rhs)
{
    const double scale = 1.f/m_dx;

    for (size_t i = 0; i < m_sizeI; i++)
    {
        for (size_t j = 0; j < m_sizeJ; j++)
        {
            if (m_materialGrid.isFluid(i,j))
            {
                double val = -scale * divergenceAt(i,j);

                val -= m_materialGrid.isSolid(i-1,j) *
                                          scale * static_cast<double>(m_fluidVelocityGrid.u(i,j) - 0);
                val += m_materialGrid.isSolid(i+1,j) *
                                          scale * static_cast<double>(m_fluidVelocityGrid.u(i+1,j) - 0);
                val -= m_materialGrid.isSolid(i,j-1) *
                                          scale * static_cast<double>(m_fluidVelocityGrid.v(i,j) - 0);
                val += m_materialGrid.isSolid(i,j+1) *
                                          scale * static_cast<double>(m_fluidVelocityGrid.v(i,j+1) - 0);

                rhs.at(linearIndex(i,j)) = val;
            }
            else
            {
                rhs.at(linearIndex(i,j)) = 0.0;
            }
        }
    }
}

void FlipSolver::calcDensityCorrectionRhs(std::vector<double> &rhs)
{
    const double scale = 1.0/m_stepDt;

    for (size_t i = 0; i < m_sizeI; i++)
    {
        for (size_t j = 0; j < m_sizeJ; j++)
        {
            if (m_materialGrid.isFluid(i,j))
            {
                rhs[linearIndex(i,j)] = scale * (1.0 - std::clamp(m_densityGrid.getAt(i,j) / m_fluidDensity,0.5,1.5));
            }
            else
            {
                rhs[linearIndex(i,j)] = 0.0;
            }
        }
    }
}

Vec3 FlipSolver::jitteredPosInCell(size_t i, size_t j)
{
    static std::uniform_real_distribution<float> dist(0.f,1.f);
    float x = (static_cast<float>(i) + dist(m_randEngine));
    float y = (static_cast<float>(j) + dist(m_randEngine));
    return Vec3(x,y);
}

void FlipSolver::countParticles()
{
    for(size_t idx = 0; idx < m_fluidVelocityGrid.linearSize(); idx++)
    {
        Index2d i2d = m_fluidVelocityGrid.index2d(idx);
        ParticleBin& bin = m_markerParticles.binForGridIdx(i2d);
        m_fluidParticleCounts.at(i2d) = 0;
        for(size_t particleIdx = 0; particleIdx < bin.size(); particleIdx++)
        {
            Vec3& position = bin.particlePosition(particleIdx);
            ssize_t i = std::floor(position.x());
            ssize_t j = std::floor(position.y());
            if(i2d.i == i && i2d.j == j)
            {
                if(m_fluidParticleCounts.at(i2d) >= 2*m_particlesPerCell)
                {
                    bin.eraseMarkerParticle(particleIdx);
                    particleIdx--;
                    continue;
                }

                m_fluidParticleCounts.at(i2d) += 1;
            }
//        if(m_materialGrid.isSolid(i,j))
//        {
//            std::cout << "Particle in solid at " << i << "," << j << '\n';
//            debug() << "Particle in solid at " << i << "," << j;
//        }
        }
    }
}

void FlipSolver::updateMaterials()
{
    for (size_t i = 0; i < m_sizeI; i++)
    {
        for (size_t j = 0; j < m_sizeJ; j++)
        {
            if(m_fluidSdf.at(i,j) < 0.f)
            {
                if(m_materialGrid.isEmpty(i,j))
                {
                    m_materialGrid.at(i,j) = FluidMaterial::FLUID;
                }
            }
            else
            {
                if(m_materialGrid.isStrictFluid(i,j))
                {
                    m_materialGrid.at(i,j) = FluidMaterial::EMPTY;
                }
            }
        }
    }
}

void FlipSolver::updateVelocityFromSolids()
{
    for (size_t i = 0; i < m_sizeI; i++)
    {
        for (size_t j = 0; j < m_sizeJ; j++)
        {
            std::vector<int> ids = validSolidNeighborIds(i,j);
            if(!ids.empty())
            {
                float avg = 0;
                for(int id : ids)
                {
                    Obstacle& obj = m_obstacles[id];
                    avg += obj.friction();
                }
                avg /= ids.size();
                if(m_materialGrid.uSampleAffectedBySolid(i,j))
                {
                    m_fluidVelocityGrid.setU(i,j,m_fluidVelocityGrid.u(i,j) * (1-avg));
                }
                if(m_materialGrid.vSampleAffectedBySolid(i,j))
                {
                    m_fluidVelocityGrid.setV(i,j,m_fluidVelocityGrid.v(i,j) * (1-avg));
                }
            }
        }
    }
}

void FlipSolver::applyPressuresToVelocityField(const std::vector<double> &pressures)
{
    std::vector<Range> ranges = ThreadPool::i()->splitRange(pressures.size());

    for(Range range : ranges)
    {
        ThreadPool::i()->enqueue(&FlipSolver::applyPressureThreadU,this,range,std::ref(pressures));
        ThreadPool::i()->enqueue(&FlipSolver::applyPressureThreadV,this,range,std::ref(pressures));
    }
    ThreadPool::i()->wait();

    // for(int i = 0; i < pressures.size(); i++)
    // {
    //     m_testGrid.data().at(i) = pressures.at(i) / 100.0;
    // }

//    if(anyNanInf(m_fluidVelocityGrid.velocityGridU().data()))
//    {
//        std::cout << "NaN or inf in U vector!\n" << std::flush;
//    }

//    if(anyNanInf(m_fluidVelocityGrid.velocityGridV().data()))
//    {
//        std::cout << "NaN or inf in V vector!\n" << std::flush;
//    }
}

void FlipSolver::applyPressureThreadU(Range range, const std::vector<double> &pressures)
{
    const double scale = m_stepDt / (m_fluidDensity * m_dx);

    for (size_t pressureIdx = range.start; pressureIdx < range.end; pressureIdx++)
    {
        Index2d i2d = index2d(pressureIdx);
        Index2d i2dim1 = Index2d(i2d.i - 1,   i2d.j);
        //Index2d i2djm1 = Index2d(i2d.m_i,       i2d.m_j - 1);
        ssize_t pressureIndexIm1 = linearIndex(i2dim1);
        double pCurrent = pressures.at(pressureIdx);
        double pIm1 = pressureIndexIm1 == -1 ? 0.0 : pressures.at(pressureIndexIm1);
        //U part
        if(m_materialGrid.isFluid(i2dim1) || m_materialGrid.isFluid(i2d))
        {
            if(m_materialGrid.isSolid(i2dim1) || m_materialGrid.isSolid(i2d))
            {
                m_fluidVelocityGrid.setU(i2d,0.f);//Solids are stationary
            }
            else
            {
                m_fluidVelocityGrid.u(i2d) -= scale * (pCurrent - pIm1);
            }
        }
        else
        {
            m_fluidVelocityGrid.setUValidity(i2d,false);
        }
    }
}

void FlipSolver::applyPressureThreadV(Range range, const std::vector<double> &pressures)
{
    const double scale = m_stepDt / (m_fluidDensity * m_dx);

    for (size_t pressureIdx = range.start; pressureIdx < range.end; pressureIdx++)
    {
        Index2d i2d = index2d(pressureIdx);
        Index2d i2djm1 = Index2d(i2d.i,   i2d.j - 1);
        //Index2d i2djm1 = Index2d(i2d.m_i,       i2d.m_j - 1);
        size_t pressureIndexJm1 = linearIndex(i2djm1);
        double pCurrent = pressures.at(pressureIdx);
        double pJm1 = pressureIndexJm1 == -1 ? 0.0 : pressures.at(pressureIndexJm1);
        //U part
        if(m_materialGrid.isFluid(i2djm1) || m_materialGrid.isFluid(i2d))
        {
            if(m_materialGrid.isSolid(i2djm1) || m_materialGrid.isSolid(i2d))
            {
                m_fluidVelocityGrid.setV(i2d,0.f);//Solids are stationary
            }
            else
            {
                m_fluidVelocityGrid.v(i2d) -= scale * (pCurrent - pJm1);
            }
        }
        else
        {
            m_fluidVelocityGrid.setVValidity(i2d,false);
        }
    }
}

Vec3 FlipSolver::rk4Integrate(Vec3 currentPosition, StaggeredVelocityGrid &grid, float dt)
{
    Vec3 k1 = dt*grid.velocityAt(currentPosition);
    Vec3 k2 = dt*grid.velocityAt(currentPosition + 0.5f*k1);
    Vec3 k3 = dt*grid.velocityAt(currentPosition + 0.5f*k2);
    Vec3 k4 = dt*grid.velocityAt(currentPosition + k3);

    return currentPosition + (1.0f/6.0f)*(k1 + 2.f*k2 + 2.f*k3 + k4);
}

void FlipSolver::gridUpdate()
{
    particleToGrid();
    m_stats.endStage(PARTICLE_TO_GRID);
    updateSdf();
    updateMaterials();
}

void FlipSolver::particleToGrid()
{
    particleVelocityToGrid();
    if(m_parameterHandlingMethod != GRID)
    {
        centeredParamsToGrid();
    }
}

void FlipSolver::applyBodyForces()
{
    const float factor = m_stepDt / m_dx;
    for (size_t i = 0; i < m_sizeI + 1; i++)
    {
        for (size_t j = 0; j < m_sizeJ + 1; j++)
        {
            if(m_fluidVelocityGrid.velocityGridU().inBounds(i,j))
            {
                m_fluidVelocityGrid.u(i,j) +=
                        factor * m_globalAcceleration.x();
            }
            if(m_fluidVelocityGrid.velocityGridV().inBounds(i,j))
            {
                m_fluidVelocityGrid.v(i,j) +=
                        factor * m_globalAcceleration.y();
            }
        }
    }
}

void FlipSolver::updateSdf()
{
    std::vector<Range> ranges = ThreadPool::i()->splitRange(m_fluidSdf.data().size());

    for(const Range& range: ranges)
    {
        ThreadPool::i()->enqueue(&FlipSolver::updateSdfThread,this,range);
    }
    ThreadPool::i()->wait();
//    float particleRadius = m_particleScale;
//    for(int i = 0; i < m_sizeI; i++)
//    {
//        for(int j = 0; j < m_sizeJ; j++)
//        {
//            float distSqrd = std::numeric_limits<float>::max();
//            Vertex centerPoint = Vertex(static_cast<float>(i) + 0.5f,
//                                        static_cast<float>(j) + 0.5f);
//            for(MarkerParticle& p : m_markerParticles)
//            {
//                float diffX = p.position.x() - centerPoint.x();
//                float diffY = p.position.y() - centerPoint.y();
//                float newDistSqrd = diffX*diffX + diffY * diffY;
//                if(newDistSqrd < distSqrd)
//                {
//                    distSqrd = newDistSqrd;
//                }
//            }

//            m_fluidSdf.at(i,j) = std::sqrt(distSqrd) - particleRadius;
//        }
//    }
}

void FlipSolver::updateSdfThread(Range range)
{
    const float particleRadius = m_particleScale;
    std::vector<float>& sdfData = m_fluidSdf.data();
//    std::vector<float>& testValues = std::get<std::vector<float>>(
//                                    m_markerParticles.getProperties(m_testValuePropertyIndex));
    for(size_t idx = range.start; idx < range.end; idx++)
    {
        Index2d i2d = m_fluidSdf.index2d(idx);
        std::array<ssize_t,9> affectingBins = m_markerParticles.binsForGridCell(i2d);
        float distSqrd = std::numeric_limits<float>::max();
        Vec3 centerPoint = Vec3(static_cast<float>(i2d.i) + 0.5f,
                                    static_cast<float>(i2d.j) + 0.5f);
        for(ssize_t binIdx : affectingBins)
        {
            if(binIdx >= 0)
            {
                ParticleBin& currentBin = m_markerParticles.binForBinIdx(binIdx);
                for(size_t particleIdx = 0; particleIdx < currentBin.size(); particleIdx++)
                {
                    //testValues[particleIdx] = binIdx >= 0;
                    Vec3 position = currentBin.particlePosition(particleIdx);
                    //m_testGrid.at(position.x(), position.y()) = 1;
                    float diffX = position.x() - centerPoint.x();
                    float diffY = position.y() - centerPoint.y();
                    float newDistSqrd = diffX*diffX + diffY * diffY;
                    if(newDistSqrd < distSqrd)
                    {
                        distSqrd = newDistSqrd;
                    }
                }
            }
        }
        sdfData[idx] = std::sqrt(distSqrd) - particleRadius;
    }
}

void FlipSolver::particleVelocityToGrid()
{
    m_fluidVelocityGrid.velocityGridU().fill(0.f);
    m_fluidVelocityGrid.velocityGridV().fill(0.f);
    m_fluidVelocityGrid.uSampleValidityGrid().fill(false);
    m_fluidVelocityGrid.vSampleValidityGrid().fill(false);

    std::vector<Range> ranges = ThreadPool::i()->splitRange(m_fluidVelocityGrid.linearSize(),1,16);

    for(Range& range : ranges)
    {
        ThreadPool::i()->enqueue(&FlipSolver::particleVelocityToGridThread,this,range);
    }
    ThreadPool::i()->wait();
}

void FlipSolver::particleVelocityToGridThread(Range r)
{
    const float threshold = 1e-9f;
    for(size_t idx = r.start; idx < r.end; idx++)
    {
        Index2d i2d = m_fluidVelocityGrid.index2d(idx);
        std::array<ssize_t,9> affectingBins = m_markerParticles.binsForGridCell(i2d);

        float uWeightTotal = 1e-10f;
        float vWeightTotal = 1e-10f;
        float uVel = 0.0;
        float vVel = 0.0;
        bool uValid = false;
        bool vValid = false;

        for(ssize_t binIdx : affectingBins)
        {
            if(binIdx >= 0)
            {
                const ParticleBin& currentBin = m_markerParticles.binForBinIdx(binIdx);
                for(size_t particleIdx = 0; particleIdx < currentBin.size(); particleIdx++)
                {
                    const Vec3 position = currentBin.particlePosition(particleIdx);
                    float weightU = simmath::quadraticBSpline(position.x() - static_cast<float>(i2d.i),
                                                              position.y() - (static_cast<float>(i2d.j) + 0.5f));

                    float weightV = simmath::quadraticBSpline(position.x() - (static_cast<float>(i2d.i) + 0.5f),
                                                              position.y() - static_cast<float>(i2d.j));
                    bool flag = weightU > threshold && weightV > threshold;
                    if(flag)
                    {
                        const Vec3 velocity = currentBin.particleVelocity(particleIdx);
                        uWeightTotal += weightU;
                        uVel += weightU * (velocity.x());
                        uValid |= true;

                        vWeightTotal += weightV;
                        vVel += weightV * (velocity.y());
                        vValid |= true;
                    }
                }
            }
        }
        m_fluidVelocityGrid.u(i2d) = uVel / uWeightTotal;
        m_fluidVelocityGrid.setUValidity(i2d,uValid);
        m_fluidVelocityGrid.v(i2d) = vVel / vWeightTotal;
        m_fluidVelocityGrid.setVValidity(i2d,vValid);
    }
}

void FlipSolver::centeredParamsToGrid()
{
    m_viscosityGrid.fill(0.f);
    m_divergenceControl.fill(0.f);
    m_knownCenteredParams.fill(false);

    std::vector<Range> ranges = ThreadPool::i()->splitRange(m_viscosityGrid.linearSize(),1,16);

    for(Range& range : ranges)
    {
        ThreadPool::i()->enqueue(&FlipSolver::centeredParamsToGridThread,this,range);
    }
    ThreadPool::i()->wait();
}

void FlipSolver::centeredParamsToGridThread(Range r)
{
    const float threshold = 1e-9f;
    for(size_t idx = r.start; idx < r.end; idx++)
    {
        Index2d i2d = m_fluidVelocityGrid.index2d(idx);
        std::array<ssize_t,9> affectingBins = m_markerParticles.binsForGridCell(i2d);

        float weightTotal = 0.0;
        float viscVal = 0.0;
        bool viscValid = false;

        for(ssize_t binIdx : affectingBins)
        {
            if(binIdx >= 0)
            {
                const ParticleBin& currentBin = m_markerParticles.binForBinIdx(binIdx);
                auto& viscVec = currentBin.particleProperties<float>(m_viscosityPropertyIndex);

                for(size_t particleIdx = 0; particleIdx < currentBin.size(); particleIdx++)
                {
                    const Vec3 position = currentBin.particlePosition(particleIdx);
                    float weightCentered = simmath::quadraticBSpline(position.x() - (i2d.i),
                                                                     position.y() - (i2d.j));
                    if(weightCentered > threshold)
                    {
                        weightTotal += weightCentered;
                        viscVal += weightCentered * viscVec.at(particleIdx);
                        viscValid = true;
                    }
                }
            }
        }

        m_viscosityGrid.at(i2d) = viscVal/weightTotal;
        m_knownCenteredParams.at(i2d) = viscValid;
    }
}

void FlipSolver::extrapolateLevelsetInside(SdfGrid &grid)
{
    ssize_t sizeI = grid.sizeI();
    ssize_t sizeJ = grid.sizeJ();
    Grid2d<int> markers(sizeI,sizeJ,std::numeric_limits<int>().max());
    std::queue<size_t> wavefront;
    for(ssize_t i = 0; i < sizeI; i++)
    {
        for(ssize_t j = 0; j < sizeJ; j++)
        {
            if(grid.at(i,j) > 0.f)
            {
                markers.at(i,j) = 0;
            }
        }
    }

    for(ssize_t i = 0; i < sizeI; i++)
    {
        for(ssize_t j = 0; j < sizeJ; j++)
        {
            if(markers.at(i,j) != 0)
            {
                for(ssize_t neighborIndex : grid.getNeighborhood(i,j))
                {
                    if(neighborIndex == -1) continue;
                    if(markers.data().at(neighborIndex) == 0)
                    {
                        markers.at(i,j) = 1;
                        wavefront.push(grid.linearIndex(i,j));
                        break;
                    }
                }
            }
        }
    }

    while(!wavefront.empty())
    {
        size_t index = wavefront.front();
        std::array<ssize_t, 8> neighbors = grid.getNeighborhood(index);
        double avg = 0;
        int count = 0;
        for(ssize_t neighborIndex : neighbors)
        {
            if(neighborIndex == -1) continue;
            if(markers.data().at(neighborIndex) < markers.data().at(index))
            {
                avg += grid.data().at(neighborIndex);
                count++;
            }
            if(markers.data().at(neighborIndex) == std::numeric_limits<int>().max() && markers.data().at(index) <= 1e6)
            {
                markers.data().at(neighborIndex) = markers.data().at(index) + 1;
                wavefront.push(neighborIndex);
            }
        }
        grid.data().at(index) = avg / count - 1.f;

        wavefront.pop();
    }
}

void FlipSolver::extrapolateLevelsetOutside(SdfGrid &grid)
{
    size_t sizeI = grid.sizeI();
    size_t sizeJ = grid.sizeJ();
    const float maxSdf = sizeI * sizeJ;
    Grid2d<int> markers(sizeI,sizeJ,std::numeric_limits<int>().max());
    std::queue<size_t> wavefront;
    for(ssize_t i = 0; i < sizeI; i++)
    {
        for(ssize_t j = 0; j < sizeJ; j++)
        {
            if(grid.at(i,j) < maxSdf)
            {
                markers.at(i,j) = 0;
            }
        }
    }

    for(ssize_t i = 0; i < sizeI; i++)
    {
        for(ssize_t j = 0; j < sizeJ; j++)
        {
            if(markers.at(i,j) != 0)
            {
                for(ssize_t neighborIndex : grid.getNeighborhood(i,j))
                {
                    if(neighborIndex == -1) continue;
                    if(markers.data().at(neighborIndex) == 0)
                    {
                        markers.at(i,j) = 1;
                        wavefront.push(grid.linearIndex(i,j));
                        break;
                    }
                }
            }
        }
    }

    while(!wavefront.empty())
    {
        size_t index = wavefront.front();
        std::array<ssize_t, 8> neighbors = grid.getNeighborhood(index);
        double avg = 0;
        int count = 0;
        for(ssize_t neighborIndex : neighbors)
        {
            if(neighborIndex == -1) continue;
            if(markers.data().at(neighborIndex) < markers.data().at(index))
            {
                avg += grid.data().at(neighborIndex);
                count++;
            }
            if(markers.data().at(neighborIndex) == std::numeric_limits<int>().max() && markers.data().at(index) <= 1e6)
            {
                markers.data().at(neighborIndex) = markers.data().at(index) + 1;
                wavefront.push(neighborIndex);
            }
        }
        grid.data().at(index) = avg / count + 1.f;

        wavefront.pop();
    }
}

float FlipSolver::maxParticleVelocity()
{
    float maxVelocitySqr = std::numeric_limits<float>::min();
    for(ParticleBin& bin: m_markerParticles.bins().data())
    {
        for(size_t particleIdx = 0; particleIdx < bin.size(); particleIdx++)
        {
            const Vec3 velocity = bin.particleVelocity(particleIdx);
            float velocitySqr = velocity.x()*velocity.x() + velocity.y()*velocity.y();
            if(velocitySqr > maxVelocitySqr) maxVelocitySqr = velocitySqr;
        }
    }

    return std::sqrt(maxVelocitySqr);
}

float FlipSolver::maxGridVelocity()
{
    float maxVelocitySqr = std::numeric_limits<float>::min();
    for(ssize_t i = 0; i < m_sizeI; i++)
    {
        for(ssize_t j = 0; j < m_sizeJ; j++)
        {
            float vU = m_fluidVelocityGrid.getU(i,j);
            float vV = m_fluidVelocityGrid.getV(i,j);
            float velocitySqr = vU*vU + vV*vV;
            if(velocitySqr > maxVelocitySqr) maxVelocitySqr = velocitySqr;
        }
    }

    return std::sqrt(maxVelocitySqr);
}

float FlipSolver::sceneScale() const
{
    return m_sceneScale;
}

float FlipSolver::lastFrameTime() const
{
    return m_frameTime;
}

float FlipSolver::avgFrameTime() const
{
    return m_avgFrameMs;
}

size_t FlipSolver::testValuePropertyIndex()
{
    return m_testValuePropertyIndex;
}

void FlipSolver::initAdditionalParameters()
{
    m_viscosityPropertyIndex = m_markerParticles.addParticleProperty<float>();
}

Vec3 FlipSolver::globalAcceleration() const
{
    return m_globalAcceleration;
}

float FlipSolver::frameDt() const
{
    return m_frameDt;
}

int FlipSolver::fps() const
{
    return m_fps;
}

int FlipSolver::maxSubsteps() const
{
    return m_maxSubsteps;
}

float FlipSolver::cflNumber() const
{
    return m_cflNumber;
}

float FlipSolver::picRatio() const
{
    return m_picRatio;
}

int FlipSolver::particlesPerCell() const
{
    return m_particlesPerCell;
}

double FlipSolver::fluidDensity() const
{
    return m_fluidDensity;
}

float FlipSolver::domainSizeJ() const
{
    return m_domainSizeJ;
}

float FlipSolver::domainSizeI() const
{
    return m_domainSizeI;
}

SimulationMethod FlipSolver::simulationMethod() const
{
    return m_simulationMethod;
}

double FlipSolver::dx() const
{
    return m_dx;
}

float FlipSolver::stepDt() const
{
    return m_stepDt;
}
