#include "flipsolver2d.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <limits>
#include <ostream>
#include <queue>
#include <thread>
#include <type_traits>

#include "Eigen/src/Core/Matrix.h"
#include "Eigen/src/Core/util/Constants.h"
#include "Eigen/src/IterativeLinearSolvers/ConjugateGradient.h"
#include "Eigen/src/IterativeLinearSolvers/IncompleteCholesky.h"
#include "Eigen/src/SparseCore/SparseMatrix.h"
#include "grid2d.h"
#include "index2d.h"
#include "linearindexable2d.h"
#include "linearsolver.h"
#include "markerparticlesystem.h"
#include "materialgrid.h"

#include "dynamicmatrix.h"
#include "logger.h"
#include "mathfuncs.h"
#include "threadpool.h"
#include "staticmatrix.h"

#include <Eigen/IterativeLinearSolvers>

FlipSolver::FlipSolver(const FlipSolverParameters *p) :
    LinearIndexable2d(p->gridSizeI, p->gridSizeJ),
    m_frameNumber(0),
    m_validVVelocitySampleCount(0),
    m_validUVelocitySampleCount(0),
    m_markerParticles(p->gridSizeI, p->gridSizeJ,5),
    m_fluidVelocityGrid(p->gridSizeI, p->gridSizeJ),
    m_savedFluidVelocityGrid(p->gridSizeI, p->gridSizeJ),
    m_materialGrid(p->gridSizeI,p->gridSizeJ, FluidMaterial::SINK),
    m_solidSdf(p->gridSizeI,p->gridSizeJ),
    m_fluidSdf(p->gridSizeI, p->gridSizeJ),
    m_knownCenteredParams(p->gridSizeI,p->gridSizeJ, false, OOBStrategy::OOB_CONST, true),
    m_viscosityGrid(p->gridSizeI,p->gridSizeJ,0.f,OOBStrategy::OOB_EXTEND, 0.f, Vertex(0.5f,0.5f)),
    m_emitterId(p->gridSizeI, p->gridSizeJ, -1),
    m_solidId(p->gridSizeI, p->gridSizeJ,-1),
    m_fluidParticleCounts(p->gridSizeI, p->gridSizeJ),
    m_divergenceControl(p->gridSizeI,p->gridSizeJ, 0.f, OOBStrategy::OOB_CONST, 0.f),
    m_densityGrid(p->gridSizeI,p->gridSizeJ, 0.f),
    m_testGrid(p->gridSizeI,p->gridSizeJ),
    m_rhs(),
    m_solverResult(),
    m_pcgSolver(m_materialGrid,10),
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
    m_simulationMethod(p->simulationMethod)
{
    Eigen::initParallel();
    Eigen::setNbThreads(ThreadPool::i()->threadCount());
    m_randEngine = std::mt19937(p->seed);
    m_testValuePropertyIndex = m_markerParticles.addParticleProperty<float>();
    m_projectTolerance = m_viscosityEnabled? 1e-6 : 1e-2;
    m_projectPreconditioner = std::make_shared<IPPreconditioner>(*dynamic_cast<LinearIndexable2d*>(this));
    m_densityPreconditioner = std::make_shared<IPPreconditioner>(*dynamic_cast<LinearIndexable2d*>(this));
    m_viscosityPreconditioner = std::make_shared<StubPreconditioner>();
    m_rhs.resize(linearSize());
    m_solverResult.resize(linearSize());
    //m_viscositySolver.setMaxIterations(200);
}

FlipSolver::~FlipSolver()
{
    ThreadPool::i()->wait();
}

void FlipSolver::project()
{
    //using precond = Eigen::IncompleteLUT<double>;
    m_rhs.resize(cellCount());
    m_solverResult.resize(cellCount());

    calcPressureRhs(m_rhs);
//    //debug() << "Calculated rhs: " << rhs;
    //Eigen::BiCGSTAB<Eigen::SparseMatrix<double>,precond> solver;

    m_solverResult = m_pressureSolver.solve(m_rhs);
    if(m_pressureSolver.info()!=Eigen::Success) {
        std::cout << "Pressure solver solving failed!\n";
        return;
    }
    std::cout << "Pressure done with " << m_pressureSolver.iterations() << " iterations\n";
//    if(!m_pcgSolver.solve(mat,m_pressures.data(),m_rhs,m_projectPreconditioner,&pd, m_pcgIterLimit, m_projectTolerance))
//    {
//        std::cout << "PCG Solver pressure: Iteration limit exhaustion!\n";
//    }
//    auto provider = getPressureMatrixElementProvider();

//    if(!m_pcgSolver.mfcgSolve(provider,m_pressures.data(),m_rhs,m_pcgIterLimit,1e-2))
//    {
//        std::cout << "PCG Solver pressure: Iteration limit exhaustion!\n";
//    }

//    if(anyNanInf(m_pressures.data()))
//    {
//        std::cout << "NaN or inf in pressures vector!\n" << std::flush;
//    }

    //debug() << "pressures = " << pressures;

    applyPressuresToVelocityField(m_solverResult);
}

LinearSolver::MatElementProvider FlipSolver::getPressureMatrixElementProvider()
{
    return std::bind(&FlipSolver::getMatFreeElementForLinIdx,this,std::placeholders::_1);
}

LinearSolver::SparseMatRowElements FlipSolver::getMatFreeElementForLinIdx(unsigned int i)
{
    double scale = m_stepDt / (m_fluidDensity * m_dx * m_dx);
    std::array<int,4> neighbors = immidiateNeighbors(static_cast<int>(i));
    std::vector<FluidMaterial>& materials = m_materialGrid.data();

    LinearSolver::SparseMatRowElements output;
    output.fill(std::pair<int, double>(0,0.0));
    output[4].first = i;

    if(fluidTest(materials[i]))
    {
        for(unsigned int i = 0; i < neighbors.size(); i++)
        {
            output[i].first = neighbors[i];
            output[i].second = -scale * fluidTest(materials[neighbors[i]]);
            output[4].second += scale * !solidTest(materials[neighbors[i]]);
        }
    }

    return output;
}

void FlipSolver::applyViscosity()
{
    //updateLinearFluidViscosityMapping();
    Eigen::VectorXd rhs;
    Eigen::VectorXd result;
    rhs.resize(m_fluidVelocityGrid.velocityGridU().linearSize() +
                 m_fluidVelocityGrid.velocityGridV().linearSize());

    result.resize(m_fluidVelocityGrid.velocityGridU().linearSize() +
                 m_fluidVelocityGrid.velocityGridV().linearSize());


    calcViscosityRhs(rhs);
    auto viscosityMatrix = getViscosityMatrix();

    m_viscositySolver.setTolerance(1e-4);
    m_viscositySolver.compute(viscosityMatrix);
    if(m_viscositySolver.info()!=Eigen::Success) {
        std::cout << "Viscosity solver decomposition failed!\n";
        return;
    }

    result = m_viscositySolver.solve(rhs);
    if(m_viscositySolver.info()!=Eigen::Success) {
        std::cout << "Viscosity solver solving failed!\n";
        return;
    }
    std::cout << "Viscosity done with " << m_viscositySolver.iterations() << " iterations\n";

    int vBaseIndex = m_fluidVelocityGrid.velocityGridU().linearSize();
    for (int i = 0; i < m_sizeI; i++)
    {
        for (int j = 0; j < m_sizeJ; j++)
        {
            int idxU = m_fluidVelocityGrid.velocityGridU().linearIndex(i,j);
            int idxV = m_fluidVelocityGrid.velocityGridV().linearIndex(i,j);
            m_fluidVelocityGrid.setU(i,j,result[idxU]);
            m_fluidVelocityGrid.setV(i,j,result[vBaseIndex + idxV]);
        }
    }

    if(anyNanInf(m_fluidVelocityGrid.velocityGridU().data()))
    {
        std::cout << "NaN or inf in U velocity after viscosity!\n" << std::flush;
    }

    if(anyNanInf(m_fluidVelocityGrid.velocityGridV().data()))
    {
        std::cout << "NaN or inf in V velocity after viscosity!\n" << std::flush;
    }
}

void FlipSolver::advect()
{
    std::vector<Range> ranges = ThreadPool::i()->splitRange(m_markerParticles.particleCount());

    for(Range& range : ranges)
    {
        ThreadPool::i()->enqueue(&FlipSolver::advectThread,this,range);
    }
    ThreadPool::i()->wait();

    //std::cout << "Advection done in max " << maxSubsteps << " substeps" << std::endl;
}

void FlipSolver::densityCorrection()
{
    updateDensityGrid();

    calcDensityCorrectionRhs(m_rhs);
    if(m_pressureSolver.info()!=Eigen::Success) {
        std::cout << "Density solver decomposition failed!\n";
        return;
    }

    m_solverResult = m_pressureSolver.solve(m_rhs);
    if(m_pressureSolver.info()!=Eigen::Success) {
        std::cout << "Density solver solving failed!\n";
        return;
    }
    std::cout << "Density done with " << m_pressureSolver.iterations() << " iterations\n";

    adjustParticlesByDensity();
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

    for(int idx = r.start; idx < r.end; idx++)
    {
        Index2d cellIdx = m_densityGrid.index2d(idx);
        for(int binIdx : m_markerParticles.binsForGridCell(cellIdx))
        {
            if(binIdx == -1)
            {
                continue;
            }

            for(size_t particleIdx : bins.data()[binIdx])
            {
                if(m_markerParticles.markedForDeath(particleIdx))
                {
                    continue;
                }

                Vertex p = m_markerParticles.particlePosition(particleIdx);
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
    std::vector<Range> ranges = ThreadPool::i()->splitRange(m_markerParticles.particleCount());
    for(Range& range : ranges)
    {
        ThreadPool::i()->enqueue(&FlipSolver::adjustParticlesByDensityThread,this,range);
    }
    ThreadPool::i()->wait();
}

void FlipSolver::adjustParticlesByDensityThread(Range r)
{
    const float scale = (m_stepDt * m_stepDt) / (m_fluidDensity * m_dx * m_dx);
    for(int idx = r.start; idx < r.end; idx++)
    {
        if(m_markerParticles.markedForDeath(idx))
        {
            continue;
        }

        Vertex& position = m_markerParticles.positions().at(idx);
        int iCorr = std::clamp(static_cast<int>(position.x() - 0.5f),0,m_sizeI);
        int jCorr = std::clamp(static_cast<int>(position.y() - 0.5f),0,m_sizeJ);
        int i = std::clamp(static_cast<int>(position.x()),0,m_sizeI);
        int j = std::clamp(static_cast<int>(position.y()),0,m_sizeJ);

        float pCurrentI = m_solverResult[linearIndex(iCorr,j)];
        float pCurrentJ = m_solverResult[linearIndex(i,jCorr)];

        float pI = m_solverResult[linearIndex(iCorr+1,j)];
        float pJ = m_solverResult[linearIndex(i,jCorr+1)];
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

void FlipSolver::advectThread(Range range)
{
    for(int i = range.start; i < range.end; i++)
    {
        Vertex& position = m_markerParticles.positions()[i];
        position = rk4Integrate(position, m_fluidVelocityGrid);
        if(m_solidSdf.interpolateAt(position.x(),position.y()) < 0.f)
        {
            position = m_solidSdf.closestSurfacePoint(position);
        }
        int pI = simmath::integr(position.x());
        int pJ = simmath::integr(position.y());
        if(!inBounds(pI,pJ) || m_materialGrid.isSink(pI,pJ))
        {
            m_markerParticles.markForDeath(i);
        }
    }
}

void FlipSolver::particleUpdate()
{
    Grid2d<float>& prevU = m_savedFluidVelocityGrid.velocityGridU();
    Grid2d<float>& prevV = m_savedFluidVelocityGrid.velocityGridV();
    for(int i = 0; i < m_markerParticles.particleCount(); i++)
    {
        Vertex &position = m_markerParticles.particlePosition(i);
        Vertex &velocity = m_markerParticles.particleVelocity(i);

        Vertex oldVelocity(prevU.interpolateAt(position.x(),position.y()),
                           prevV.interpolateAt(position.x(),position.y()));
        Vertex newVelocity = m_fluidVelocityGrid.velocityAt(position) ;
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

void FlipSolver::afterTransfer()
{
    for (int i = 0; i < m_sizeI; i++)
    {
        for (int j = 0; j < m_sizeJ; j++)
        {
            if(m_materialGrid.isSource(i,j))
            {
                int emitterId = m_emitterId.at(i,j);
                m_viscosityGrid.at(i,j) = m_sources[emitterId].viscosity();
//                if(m_sources[emitterId].velocityTransfer())
//                {
                    m_fluidVelocityGrid.setU(i,j,m_sources[emitterId].velocity().x() / m_dx);
                    m_fluidVelocityGrid.setV(i,j,m_sources[emitterId].velocity().y() / m_dx);
                    m_fluidVelocityGrid.setUValidity(i,j,true);
                    m_fluidVelocityGrid.setVValidity(i,j,true);
//                }
            }
        }
    }
}

void FlipSolver::step()
{
    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;
    using std::chrono::milliseconds;
    static double sum = 0.0;
    static double count = 0.0;
    auto t1 = high_resolution_clock::now();
    advect();
    auto t2 = high_resolution_clock::now();
    m_pressureMatrix = getPressureProjectionMatrix();
    m_pressureSolver.compute(m_pressureMatrix);
    if(m_pressureSolver.info()!=Eigen::Success) {
        std::cout << "Pressure solver decomposition failed!\n";
        return;
    }

    m_pressureSolver.setTolerance(1e-4);
    densityCorrection();
    m_markerParticles.pruneParticles();
    m_markerParticles.rebinParticles();
    particleToGrid();

    updateSdf();

    //updateLinearFluidViscosityMapping();
    updateMaterials();
    afterTransfer();
    extrapolateLevelsetInside(m_fluidSdf);
    m_fluidVelocityGrid.extrapolate(10);

    m_savedFluidVelocityGrid = m_fluidVelocityGrid;
    applyBodyForces();

    m_pressureSolver.setTolerance(1e-4);
    project();

    updateVelocityFromSolids();
    if(m_viscosityEnabled)
    {
        applyViscosity();
        project();
    }
    m_fluidVelocityGrid.extrapolate(10);
    particleUpdate();
    countParticles();
    reseedParticles();

    duration<double, std::milli> ms_double = t2 - t1;
    sum += ms_double.count();
    count+=1.0;
    //m_frameTime = sum / count;
    std::cout << ms_double.count() << "ms\n";
}

void FlipSolver::stepFrame()
{
    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;
    using std::chrono::milliseconds;
    static double msSum = 0.0;
    auto t1 = high_resolution_clock::now();
    if(m_frameNumber == 0)
    {
        firstFrameInit();
    }
    m_testGrid.fill(0.f);
    float substepTime = 0.f;
    bool finished = false;
    int substepCount = 0;

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
        substepTime += maxSubstepSize;
        substepCount++;
        if(substepCount > 50) break;
    }
    auto t2 = high_resolution_clock::now();
    duration<double, std::milli> ms = t2 - t1;
    //std::cout << "Frame took" << ms.count() << "ms\n";
    m_frameTime = ms.count();
    msSum += ms.count();
    std::cout << "Frame done in " << substepCount << " substeps" << std::endl;
    m_frameNumber++;
    m_avgFrameMs = msSum / m_frameNumber;
}

void FlipSolver::updateSolids()
{
    for (int i = 0; i < m_sizeI; i++)
    {
        for (int j = 0; j < m_sizeJ; j++)
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
    float dx = m_dx;
    for (int i = 0; i < m_sizeI; i++)
    {
        for (int j = 0; j < m_sizeJ; j++)
        {
            for(int emitterIdx = 0; emitterIdx < m_sources.size(); emitterIdx++)
            {
                Emitter& e = m_sources[emitterIdx];
                if(e.geometry().signedDistance(i*dx,j*dx) <= 0.f)
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
    for (int i = 0; i < m_sizeI; i++)
    {
        for (int j = 0; j < m_sizeJ; j++)
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
    for (int i = 0; i < m_sizeI; i++)
    {
        for (int j = 0; j < m_sizeJ; j++)
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

int FlipSolver::gridSizeI()
{
    return m_sizeI;
}

int FlipSolver::gridSizeJ()
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
    for (int pIndex = 0; pIndex < m_markerParticles.particleCount(); pIndex++)
    {
        Vertex pos = m_markerParticles.particlePosition(pIndex);
        int i = pos.x();
        int j = pos.y();

        if(m_fluidParticleCounts.at(i,j) > 2*m_particlesPerCell)
        {
            m_markerParticles.markForDeath(pIndex);
            m_fluidParticleCounts.at(i,j) -= 1;
            //pIndex--;
            continue;
        }
    }

    std::vector<float>& particleViscosities = m_markerParticles.particleProperties<float>
                                              (m_viscosityPropertyIndex);

    for (int i = 0; i < m_sizeI; i++)
    {
        for (int j = 0; j < m_sizeJ; j++)
        {
            int particleCount = m_fluidParticleCounts.at(i,j);
            if(particleCount > 20)
            {
                std::cout << "too many particles " << particleCount << " at " << i << ' ' << j;
            }
            //std::cout << particleCount << " at " << i << " , " << j << std::endl;
            int additionalParticles = (m_particlesPerCell / 2) - particleCount;
            if(additionalParticles <= 0)
            {
                continue;
            }
            if(m_materialGrid.isSource(i,j))
            {
                for(int p = 0; p < additionalParticles; p++)
                {
                    Vertex pos = jitteredPosInCell(i,j);
                    Vertex velocity = m_fluidVelocityGrid.velocityAt(pos);
                    //Vertex velocity = Vertex();
                    int emitterId = m_emitterId.at(i,j);
                    float viscosity = m_sources[emitterId].viscosity();
                    size_t pIdx = m_markerParticles.addMarkerParticle(pos,velocity);
                    particleViscosities[pIdx] = viscosity;
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
    std::vector<float>& particleViscosities = m_markerParticles.particleProperties<float>
                                              (m_viscosityPropertyIndex);
    for (int i = 0; i < m_sizeI; i++)
    {
        for (int j = 0; j < m_sizeJ; j++)
        {
            if(m_materialGrid.isStrictFluid(i,j))
            {
                for(int p = 0; p < m_particlesPerCell; p++)
                {
                    Vertex pos = jitteredPosInCell(i,j);
                    Vertex velocity = m_fluidVelocityGrid.velocityAt(pos);
                    float viscosity = m_viscosityGrid.interpolateAt(pos);
                    size_t pIdx = m_markerParticles.addMarkerParticle(pos,velocity);
                    particleViscosities.at(pIdx) = viscosity;
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

double FlipSolver::divergenceAt(int i, int j)
{
    return static_cast<double>(m_fluidVelocityGrid.u(i+1,j) - m_fluidVelocityGrid.u(i,j)
                               + m_fluidVelocityGrid.v(i,j+1) - m_fluidVelocityGrid.v(i,j)
                               + m_divergenceControl.at(i,j));
}

std::vector<int> FlipSolver::validSolidNeighborIds(int i, int j)
{
    std::vector<int> output;
    output.reserve(8);
    for (int iOffset = -1; iOffset < 2; iOffset++)
    {
        for (int jOffset = -1; jOffset < 2; jOffset++)
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

Eigen::SparseMatrix<double,Eigen::RowMajor> FlipSolver::getPressureProjectionMatrix()
{
    Eigen::SparseMatrix<double,Eigen::RowMajor> output = Eigen::SparseMatrix<double>();
    output.resize(cellCount(),cellCount());
    output.reserve(Eigen::VectorXi::Constant(cellCount(),6));

    const double scale = m_stepDt / (m_fluidDensity * m_dx * m_dx);

    LinearIndexable2d& indexer = *dynamic_cast<LinearIndexable2d*>(this);

    for(int i = 0; i < m_sizeI; i++)
    {
        for(int j = 0; j < m_sizeJ; j++)
        {
            const int linIdx = indexer.linearIndex(i,j);
            if(!m_materialGrid.isFluid(i,j))
            {
                output.insert(linIdx,linIdx) = 1.f;
                continue;
            }

            const int linIdxAx = indexer.linearIdxOfOffset(linIdx,1,0);
            const int linIdxAy = indexer.linearIdxOfOffset(linIdx,0,1);

            double diag = 0.0;
            //X Neighbors
            if(m_materialGrid.isFluid(i-1,j))
            {
                diag += scale;
            }else if(m_materialGrid.isEmpty(i-1,j))
            {
                diag += scale;
            }

            if(m_materialGrid.isFluid(i+1,j))
            {
                diag += scale;
                if(inBounds(linIdxAx))
                {
                    output.insert(linIdxAx,linIdx) = -scale;
                    output.insert(linIdx,linIdxAx) = -scale;
                }
            } else if(m_materialGrid.isEmpty(i+1,j))
            {
                diag += scale;
            }

            //Y Neighbors
            if(m_materialGrid.isFluid(i,j-1))
            {
                diag += scale;
            }else if(m_materialGrid.isEmpty(i,j-1))
            {
                diag += scale;
            }

            if(m_materialGrid.isFluid(i,j+1))
            {
                diag += scale;
                if(inBounds(linIdxAy))
                {
                    output.insert(linIdx,linIdxAy) = -scale;
                    output.insert(linIdxAy,linIdx) = -scale;
                }
            } else if(m_materialGrid.isEmpty(i,j+1))
            {
                diag += scale;
            }

            output.insert(linIdx,linIdx) = diag;
        }
    }

//    for(int i = 0; i <  m_sizeI; i++)
//    {
//        for(int j = 0; j <  m_sizeJ; j++)
//        {
//            if(m_materialGrid.isFluid(i,j))
//            {
//                //X Neighbors
//                if(m_materialGrid.isFluid(i-1,j))
//                {
//                    output.addToAdiag(i,j,scale,indexer);
//                }else if(m_materialGrid.isEmpty(i-1,j))
//                {
//                    output.addToAdiag(i,j,scale,  indexer);
//                }

//                if(m_materialGrid.isFluid(i+1,j))
//                {
//                    output.addToAdiag(i,j,scale,  indexer);
//                    output.setAx(i,j,-scale,  indexer);
//                } else if(m_materialGrid.isEmpty(i+1,j))
//                {
//                    output.addToAdiag(i,j,scale,  indexer);
//                }

//                //Y Neighbors
//                if(m_materialGrid.isFluid(i,j-1))
//                {
//                    output.addToAdiag(i,j,scale,  indexer);
//                }else if(m_materialGrid.isEmpty(i,j-1))
//                {
//                    output.addToAdiag(i,j,scale,  indexer);
//                }

//                if(m_materialGrid.isFluid(i,j+1))
//                {
//                    output.addToAdiag(i,j,scale,  indexer);
//                    output.setAy(i,j,-scale,  indexer);
//                } else if(m_materialGrid.isEmpty(i,j+1))
//                {
//                    output.addToAdiag(i,j,scale,  indexer);
//                }
//            }
//        }
//    }

    output.makeCompressed();

    return output;
}

Eigen::SparseMatrix<double,Eigen::RowMajor> FlipSolver::getViscosityMatrix()
{
    const size_t size = m_fluidVelocityGrid.velocityGridU().linearSize() +
                  m_fluidVelocityGrid.velocityGridV().linearSize();

    Eigen::SparseMatrix<double,Eigen::RowMajor> output = Eigen::SparseMatrix<double>();
    output.resize(size,size);
    output.reserve(Eigen::VectorXi::Constant(size,10));

    const float scaleTwoDt = 2*m_stepDt / (m_dx * m_dx);
    const float scaleTwoDx = m_stepDt / (2 * m_dx * m_dx);
    const int vBaseIndex = m_fluidVelocityGrid.velocityGridU().linearSize();

    LinearIndexable2d& uIndexer = m_fluidVelocityGrid.velocityGridU();
    LinearIndexable2d& vIndexer = m_fluidVelocityGrid.velocityGridV();

    for(int i = 0; i < m_sizeI+1; i++)
    {
        for(int j = 0; j < m_sizeJ+1; j++)
        {   
            int idxU = m_fluidVelocityGrid.velocityGridU().linearIndex(i,j);
            int idxV = m_fluidVelocityGrid.velocityGridV().linearIndex(i,j);

            if(idxU != -1)
            {
                float fi = static_cast<float>(i);
                float fj = static_cast<float>(j);

                //U component
                output.coeffRef(idxU,idxU) += m_fluidDensity;

                int uImOneLinearIdx = uIndexer.linearIndex(i-1,j);

                if(uImOneLinearIdx != -1)
                {
                    output.coeffRef(idxU,
                                 uImOneLinearIdx) += -scaleTwoDt * m_viscosityGrid.getAt(i-1,j);

                    output.coeffRef(idxU,
                                 idxU) += scaleTwoDt * m_viscosityGrid.getAt(i-1,j);
                }

                int uIpOneLinearIdx = uIndexer.linearIndex(i+1,j);

                if(uIpOneLinearIdx != -1)
                {
                    output.coeffRef(idxU,
                                 uIpOneLinearIdx) += -scaleTwoDt * m_viscosityGrid.getAt(i,j);

                    output.coeffRef(idxU,
                                 idxU) += scaleTwoDt * m_viscosityGrid.getAt(i,j);
                }

                int uJmOneLinearIdx = uIndexer.linearIndex(i,j-1);
                int vImOneLinearIdx = vIndexer.linearIndex(i-1,j);

                if(uJmOneLinearIdx != -1
                        && idxV != -1
                        && vImOneLinearIdx != -1)
                {
                    float lerpedViscosity = m_viscosityGrid.interpolateAt(fi-0.5f,fj-0.5f);
                    //lerpedViscosity = tempVisc;
                    output.coeffRef(idxU,
                                 uJmOneLinearIdx) += -scaleTwoDx * lerpedViscosity;

                    output.coeffRef(idxU,
                                    vBaseIndex + idxV) += scaleTwoDx * lerpedViscosity;

                    output.coeffRef(idxU,
                                    vBaseIndex + vImOneLinearIdx) += -scaleTwoDx * lerpedViscosity;

                    output.coeffRef(idxU,idxU) += scaleTwoDx * lerpedViscosity;
                }

                int uJpOneLinearIdx = uIndexer.linearIndex(i,j+1);
                int vJpOneLinearIdx = vIndexer.linearIndex(i,j+1);
                int vImOneJpOneLinearIdx = vIndexer.linearIndex(i-1,j+1);

                if(uJpOneLinearIdx != -1
                        && vJpOneLinearIdx != -1
                        && vImOneJpOneLinearIdx != -1)
                {
                    float lerpedViscosity = m_viscosityGrid.interpolateAt(fi-0.5f,fj+0.5f);
                    //lerpedViscosity = tempVisc;
                    output.coeffRef(idxU,
                                    uJpOneLinearIdx) += -scaleTwoDx * lerpedViscosity;

                    output.coeffRef(idxU,
                                    vBaseIndex + vJpOneLinearIdx) += -scaleTwoDx * lerpedViscosity;

                    output.coeffRef(idxU,
                                    vBaseIndex + vImOneJpOneLinearIdx) += scaleTwoDx * lerpedViscosity;

                    output.coeffRef(idxU,idxU) += scaleTwoDx * lerpedViscosity;
                }
            }

            //V component
            if(idxV != -1)
            {
                float fi = static_cast<float>(i);
                float fj = static_cast<float>(j);
                int vBaseIndex = uIndexer.linearSize();
                output.coeffRef(vBaseIndex + idxV,vBaseIndex + idxV) += m_fluidDensity;

                int vJmOneLinearIdx = vIndexer.linearIndex(i,j-1);

                if(vJmOneLinearIdx != -1)
                {
                    output.coeffRef(vBaseIndex + idxV,
                                    vBaseIndex + vJmOneLinearIdx) += -scaleTwoDt * m_viscosityGrid.getAt(i,j-1);

                    output.coeffRef(vBaseIndex + idxV,
                                 vBaseIndex + idxV) += scaleTwoDt * m_viscosityGrid.getAt(i,j-1);
                }

                int vJpOneLinearIdx = vIndexer.linearIndex(i,j+1);

                if(vJpOneLinearIdx != -1)
                {
                    output.coeffRef(vBaseIndex + idxV,
                                    vBaseIndex + vJpOneLinearIdx)
                        += -scaleTwoDt * m_viscosityGrid.getAt(i,j);

                    output.coeffRef(vBaseIndex + idxV,
                                    vBaseIndex + idxV)
                        += scaleTwoDt * m_viscosityGrid.getAt(i,j);
                }

                int uJmOneLinearIdx = uIndexer.linearIndex(i,j-1);
                int vImOneLinearIdx = vIndexer.linearIndex(i-1,j);

                if(idxU != -1
                        && uJmOneLinearIdx != -1
                        && vImOneLinearIdx != -1)
                {
                    float lerpedViscosity = m_viscosityGrid.interpolateAt(fi-0.5f,fj-0.5f);
                    //lerpedViscosity = tempVisc;

                    output.coeffRef(vBaseIndex + idxV,
                                 idxU) +=
                                 scaleTwoDx * lerpedViscosity;

                    output.coeffRef(vBaseIndex + idxV,
                                 uJmOneLinearIdx) +=
                                 -scaleTwoDx * lerpedViscosity;

                    output.coeffRef(vBaseIndex + idxV,
                                 vBaseIndex + vImOneLinearIdx) +=
                                 -scaleTwoDx * lerpedViscosity;

                    output.coeffRef(vBaseIndex + idxV,
                                 vBaseIndex + idxV) +=
                                 scaleTwoDx * lerpedViscosity;
                }

                int uIpOneLinearIdx = uIndexer.linearIndex(i+1,j);
                int uIpOneJmOneLinearIdx = uIndexer.linearIndex(i+1,j-1);
                int vIpOneLinearIdx = vIndexer.linearIndex(i+1,j);

                if(uIpOneLinearIdx != -1
                        && uIpOneJmOneLinearIdx != -1
                        && vIpOneLinearIdx != -1)
                {
                    float lerpedViscosity = m_viscosityGrid.interpolateAt(fi+0.5f,fj-0.5f);
                    //lerpedViscosity = tempVisc;

                    output.coeffRef(vBaseIndex + idxV,
                                 uIpOneLinearIdx) +=
                        -scaleTwoDx * lerpedViscosity;

                    output.coeffRef(vBaseIndex + idxV,
                                    uIpOneJmOneLinearIdx) +=
                                 scaleTwoDx * lerpedViscosity;

                    output.coeffRef(vBaseIndex + idxV,
                                 vBaseIndex + vIpOneLinearIdx) +=
                                 -scaleTwoDx * lerpedViscosity;

                    output.coeffRef(vBaseIndex + idxV,
                                    vBaseIndex + idxV) +=
                                 scaleTwoDx * lerpedViscosity;
                }
            }
        }
    }

    output.makeCompressed();

    return output;
}

size_t FlipSolver::particleCount()
{
    return m_markerParticles.particleCount();
}

int FlipSolver::cellCount()
{
    return m_sizeI * m_sizeJ;
}

void FlipSolver::calcPressureRhs(Eigen::VectorXd &rhs)
{
    double scale = 1.f/m_dx;

    for (int i = 0; i < m_sizeI; i++)
    {
        for (int j = 0; j < m_sizeJ; j++)
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

                rhs.coeffRef(linearIndex(i,j)) = val;
            }
            else
            {
                rhs.coeffRef(linearIndex(i,j)) = 0.0;
            }
        }
    }
}

void FlipSolver::calcViscosityRhs(Eigen::VectorXd &rhs)
{
    LinearIndexable2d& uIndexer = m_fluidVelocityGrid.velocityGridU();
    LinearIndexable2d& vIndexer = m_fluidVelocityGrid.velocityGridV();
    int vBaseIndex = uIndexer.linearSize();

    for (int i = 0; i < m_sizeI+1; i++)
    {
        for (int j = 0; j < m_sizeJ+1; j++)
        {
            int idxU = uIndexer.linearIndex(i,j);
            int idxV = vIndexer.linearIndex(i,j);

            if(idxU != -1)
            {
                float u = m_fluidVelocityGrid.getU(i,j);
                rhs[idxU] = m_fluidDensity * u;
            }

            if(idxV != -1)
            {
                float v = m_fluidVelocityGrid.getV(i,j);
                rhs[vBaseIndex + idxV] = m_fluidDensity * v;
            }
        }
    }
}

void FlipSolver::calcDensityCorrectionRhs(Eigen::VectorXd &rhs)
{
    const double scale = 1.0/m_stepDt;

    for (int i = 0; i < m_sizeI; i++)
    {
        for (int j = 0; j < m_sizeJ; j++)
        {
            if (m_materialGrid.isFluid(i,j))
            {
                rhs[linearIndex(i,j)] = scale * (1.0 - std::clamp(m_densityGrid.getAt(i,j) / m_fluidDensity,0.5,1.5));
            }
        }
    }
}

Vertex FlipSolver::jitteredPosInCell(int i, int j)
{
    static std::uniform_real_distribution<float> dist(0.f,1.f);
    float x = (static_cast<float>(i) + dist(m_randEngine));
    float y = (static_cast<float>(j) + dist(m_randEngine));
    return Vertex(x,y);
}

void FlipSolver::countParticles()
{
    for(size_t idx = 0; idx < m_fluidVelocityGrid.linearSize(); idx++)
    {
        Index2d i2d = m_fluidVelocityGrid.index2d(idx);
        MarkerParticleSystem::ParticleBin& bin = m_markerParticles.binForGridIdx(i2d);
        m_fluidParticleCounts.at(i2d) = 0;
        for(size_t particleIdx : bin)
        {
            Vertex& position = m_markerParticles.positions()[particleIdx];
            int i = std::floor(position.x());
            int j = std::floor(position.y());
            if(i2d.i == i && i2d.j == j)
            {
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
    for (int i = 0; i < m_sizeI; i++)
    {
        for (int j = 0; j < m_sizeJ; j++)
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
    for (int i = 0; i < m_sizeI; i++)
    {
        for (int j = 0; j < m_sizeJ; j++)
        {
            std::vector<int> ids = validSolidNeighborIds(i,j);
            if(!ids.empty())
            {
                float avg = 0;
                for(int& id : ids)
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

void FlipSolver::applyPressuresToVelocityField(Eigen::VectorXd &pressures)
{
    std::vector<Range> ranges = ThreadPool::i()->splitRange(pressures.size());

    for(Range range : ranges)
    {
        ThreadPool::i()->enqueue(&FlipSolver::applyPressureThreadU,this,range,pressures);
        ThreadPool::i()->enqueue(&FlipSolver::applyPressureThreadV,this,range,pressures);
    }
    ThreadPool::i()->wait();

//    if(anyNanInf(m_fluidVelocityGrid.velocityGridU().data()))
//    {
//        std::cout << "NaN or inf in U vector!\n" << std::flush;
//    }

//    if(anyNanInf(m_fluidVelocityGrid.velocityGridV().data()))
//    {
//        std::cout << "NaN or inf in V vector!\n" << std::flush;
//    }
}

void FlipSolver::applyPressureThreadU(Range range, const Eigen::VectorXd &pressures)
{
    const double scale = m_stepDt / (m_fluidDensity * m_dx);

    for (unsigned int pressureIdx = range.start; pressureIdx < range.end; pressureIdx++)
    {
        Index2d i2d = index2d(pressureIdx);
        Index2d i2dim1 = Index2d(i2d.i - 1,   i2d.j);
        //Index2d i2djm1 = Index2d(i2d.m_i,       i2d.m_j - 1);
        int pressureIndexIm1 = linearIndex(i2dim1);
        double pCurrent = pressures.coeffRef(pressureIdx);
        double pIm1 = pressureIndexIm1 == -1 ? 0.0 : pressures.coeffRef(pressureIndexIm1);
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

void FlipSolver::applyPressureThreadV(Range range, const Eigen::VectorXd &pressures)
{
    const double scale = m_stepDt / (m_fluidDensity * m_dx);

    for (unsigned int pressureIdx = range.start; pressureIdx < range.end; pressureIdx++)
    {
        Index2d i2d = index2d(pressureIdx);
        Index2d i2djm1 = Index2d(i2d.i,   i2d.j - 1);
        //Index2d i2djm1 = Index2d(i2d.m_i,       i2d.m_j - 1);
        int pressureIndexJm1 = linearIndex(i2djm1);
        double pCurrent = pressures.coeffRef(pressureIdx);
        double pJm1 = pressureIndexJm1 == -1 ? 0.0 : pressures.coeffRef(pressureIndexJm1);
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

Vertex FlipSolver::rk4Integrate(Vertex currentPosition, StaggeredVelocityGrid &grid)
{
    float factor = m_stepDt;
    Vertex k1 = factor*grid.velocityAt(currentPosition);
    Vertex k2 = factor*grid.velocityAt(currentPosition + 0.5f*k1);
    Vertex k3 = factor*grid.velocityAt(currentPosition + 0.5f*k2);
    Vertex k4 = factor*grid.velocityAt(currentPosition + k3);

    return currentPosition + (1.0f/6.0f)*(k1 + 2.f*k2 + 2.f*k3 + k4);
}

void FlipSolver::particleToGrid()
{
    particleVelocityToGrid();
    centeredParamsToGrid();
}

void FlipSolver::applyBodyForces()
{
    const float factor = m_stepDt / m_dx;
    for (int i = 0; i < m_sizeI + 1; i++)
    {
        for (int j = 0; j < m_sizeJ + 1; j++)
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
        std::array<int,9> affectingBins = m_markerParticles.binsForGridCell(i2d);
        float distSqrd = std::numeric_limits<float>::max();
        Vertex centerPoint = Vertex(static_cast<float>(i2d.i) + 0.5f,
                                    static_cast<float>(i2d.j) + 0.5f);
        for(int binIdx : affectingBins)
        {
            if(binIdx >= 0)
            {
                for(size_t particleIdx : m_markerParticles.binForBinIdx(binIdx))
                {
                    //testValues[particleIdx] = binIdx >= 0;
                    Vertex position = m_markerParticles.particlePosition(particleIdx);
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
    Grid2d<float> uWeights(m_fluidVelocityGrid.velocityGridU().sizeI(),
                           m_fluidVelocityGrid.velocityGridU().sizeJ(),1e-10f);
    Grid2d<float> vWeights(m_fluidVelocityGrid.velocityGridV().sizeI(),
                           m_fluidVelocityGrid.velocityGridV().sizeJ(),1e-10f);
    m_fluidVelocityGrid.velocityGridU().fill(0.f);
    m_fluidVelocityGrid.velocityGridV().fill(0.f);
    m_fluidVelocityGrid.uSampleValidityGrid().fill(false);
    m_fluidVelocityGrid.vSampleValidityGrid().fill(false);

    std::vector<Range> ranges = ThreadPool::i()->splitRange(m_fluidVelocityGrid.linearSize(),1,8);

    for(Range& range : ranges)
    {
        ThreadPool::i()->enqueue(&FlipSolver::particleVelocityToGridThread,this,range,
                                                                    std::ref(uWeights), std::ref(vWeights));
    }
    ThreadPool::i()->wait();
}

void FlipSolver::particleVelocityToGridThread(Range r, Grid2d<float> &uWeights, Grid2d<float> &vWeights)
{
    for(size_t idx = r.start; idx < r.end; idx++)
    {
        Index2d i2d = m_fluidVelocityGrid.index2d(idx);
        std::array<int,9> affectingBins = m_markerParticles.binsForGridCell(i2d);
        for(int binIdx : affectingBins)
        {
            if(binIdx >= 0)
            {
                for(size_t particleIdx : m_markerParticles.binForBinIdx(binIdx))
                {
                    Vertex position = m_markerParticles.particlePosition(particleIdx);
                    float weightU = simmath::quadraticBSpline(position.x() - static_cast<float>(i2d.i),
                                                              position.y() - (static_cast<float>(i2d.j) + 0.5f));

                    float weightV = simmath::quadraticBSpline(position.x() - (static_cast<float>(i2d.i) + 0.5f),
                                                              position.y() - static_cast<float>(i2d.j));
                    if(std::abs(weightU) > 1e-9f && std::abs(weightV) > 1e-9f)
                    {
                        Vertex velocity = m_markerParticles.particleVelocity(particleIdx);
                        uWeights.at(i2d.i,i2d.j) += weightU;
                        m_fluidVelocityGrid.u(i2d.i,i2d.j) += weightU * (velocity.x());
                        m_fluidVelocityGrid.setUValidity(i2d.i,i2d.j,true);

                        vWeights.at(i2d.i,i2d.j) += weightV;
                        m_fluidVelocityGrid.v(i2d.i,i2d.j) += weightV * (velocity.y());
                        m_fluidVelocityGrid.setVValidity(i2d.i,i2d.j,true);
                    }
                }
            }
        }
        if(m_fluidVelocityGrid.velocityGridU().inBounds(i2d))
        {
            m_fluidVelocityGrid.u(i2d) /= uWeights.at(i2d);
        }
        if(m_fluidVelocityGrid.velocityGridV().inBounds(i2d))
        {
            m_fluidVelocityGrid.v(i2d) /= vWeights.at(i2d);
        }
    }
}

void FlipSolver::centeredParamsToGrid()
{
    Grid2d<float> centeredWeights(m_sizeI,m_sizeJ,1e-10f);

    m_viscosityGrid.fill(0.f);
    m_divergenceControl.fill(0.f);
    m_knownCenteredParams.fill(false);

    std::vector<Range> ranges = ThreadPool::i()->splitRange(m_viscosityGrid.linearSize(),1,8);

    for(Range& range : ranges)
    {
        ThreadPool::i()->enqueue(&FlipSolver::centeredParamsToGridThread,this,range,
                                 std::ref(centeredWeights));
    }
    ThreadPool::i()->wait();
}

void FlipSolver::centeredParamsToGridThread(Range r, Grid2d<float> &cWeights)
{
    std::vector<float>& particleViscosities = m_markerParticles.particleProperties<float>
                                              (m_viscosityPropertyIndex);

    for(size_t idx = r.start; idx < r.end; idx++)
    {
        Index2d i2d = m_fluidVelocityGrid.index2d(idx);
        std::array<int,9> affectingBins = m_markerParticles.binsForGridCell(i2d);
        for(int binIdx : affectingBins)
        {
            if(binIdx >= 0)
            {
                for(size_t particleIdx : m_markerParticles.binForBinIdx(binIdx))
                {
                    Vertex position = m_markerParticles.particlePosition(particleIdx);
                    float weightCentered = simmath::quadraticBSpline(position.x() - (i2d.i),
                                                                     position.y() - (i2d.j));
                    if(std::abs(weightCentered) > 1e-6f)
                    {
                        cWeights.at(i2d.i,i2d.j) += weightCentered;
                        m_viscosityGrid.at(i2d.i,i2d.j) += weightCentered * particleViscosities.at(particleIdx);
                        m_knownCenteredParams.at(i2d.i,i2d.j) = true;
                    }
                }
            }
        }

        if(m_knownCenteredParams.inBounds(i2d))
        {
            if(m_knownCenteredParams.at(i2d))
            {
                m_viscosityGrid.at(i2d) /= cWeights.at(i2d);
                //m_testGrid.at(i2d) = m_viscosityGrid.at(i2d);
            }
            //m_testGrid.at(i2d) = m_viscosityGrid.at(i2d);
        }
    }
}

void FlipSolver::extrapolateLevelsetInside(SdfGrid &grid)
{
    int sizeI = grid.sizeI();
    int sizeJ = grid.sizeJ();
    Grid2d<int> markers(sizeI,sizeJ,std::numeric_limits<int>().max());
    std::queue<int> wavefront;
    for(int i = 0; i < sizeI; i++)
    {
        for(int j = 0; j < sizeJ; j++)
        {
            if(grid.at(i,j) > 0.f)
            {
                markers.at(i,j) = 0;
            }
        }
    }

    for(int i = 0; i < sizeI; i++)
    {
        for(int j = 0; j < sizeJ; j++)
        {
            if(markers.at(i,j) != 0)
            {
                for(int neighborIndex : grid.getNeighborhood(i,j))
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
        int index = wavefront.front();
        std::array<int, 8> neighbors = grid.getNeighborhood(index);
        double avg = 0;
        int count = 0;
        for(int neighborIndex : neighbors)
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
    int sizeI = grid.sizeI();
    int sizeJ = grid.sizeJ();
    const float maxSdf = sizeI * sizeJ;
    Grid2d<int> markers(sizeI,sizeJ,std::numeric_limits<int>().max());
    std::queue<int> wavefront;
    for(int i = 0; i < sizeI; i++)
    {
        for(int j = 0; j < sizeJ; j++)
        {
            if(grid.at(i,j) < maxSdf)
            {
                markers.at(i,j) = 0;
            }
        }
    }

    for(int i = 0; i < sizeI; i++)
    {
        for(int j = 0; j < sizeJ; j++)
        {
            if(markers.at(i,j) != 0)
            {
                for(int neighborIndex : grid.getNeighborhood(i,j))
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
        int index = wavefront.front();
        std::array<int, 8> neighbors = grid.getNeighborhood(index);
        double avg = 0;
        int count = 0;
        for(int neighborIndex : neighbors)
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
    for(const Vertex& velocity : m_markerParticles.velocities())
    {
        float velocitySqr = velocity.x()*velocity.x() + velocity.y()*velocity.y();
        if(velocitySqr > maxVelocitySqr) maxVelocitySqr = velocitySqr;
    }

    return std::sqrt(maxVelocitySqr);
}

float FlipSolver::maxGridVelocity()
{
    float maxVelocitySqr = std::numeric_limits<float>::min();
    for(int i = 0; i < m_sizeI; i++)
    {
        for(int j = 0; j < m_sizeJ; j++)
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

Vertex FlipSolver::globalAcceleration() const
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
