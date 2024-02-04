#include "nbflipsolver.h"
#include "flipsolver2d.h"
#include "sdfgrid.h"

#include "staggeredvelocitygrid.h"
#include <algorithm>
#include <iostream>

NBFlipSolver::NBFlipSolver(const NBFlipParameters *p):
    FlipSolver(p),
    m_advectedVelocity(p->gridSizeI, p->gridSizeJ),
    m_advectedSdf(p->gridSizeI, p->gridSizeJ),
    m_uWeights(m_fluidVelocityGrid.velocityGridU().sizeI(),
            m_fluidVelocityGrid.velocityGridU().sizeJ(),1e-10f),
    m_vWeights(m_fluidVelocityGrid.velocityGridV().sizeI(),
            m_fluidVelocityGrid.velocityGridV().sizeJ(),1e-10f),
    m_narrowBand(-3.f),
    m_combinationBand(-2.f),
    m_resamplingBand(-1.f)
{
    m_projectTolerance = 1e-6;
    m_pressureSolver.setTolerance(m_projectTolerance);
}

void NBFlipSolver::step()
{
    advect();
    m_stats.endStage(ADVECTION);
    //m_testGrid = m_viscosityGrid;

    m_markerParticles.pruneParticles();
    m_markerParticles.rebinParticles();
    m_stats.endStage(PARTICLE_REBIN);
    particleToGrid();
    m_fluidVelocityGrid.extrapolate(10);
    m_stats.endStage(PARTICLE_TO_GRID);
    m_savedFluidVelocityGrid = m_fluidVelocityGrid;

    updateSdf();
    //updateLinearFluidViscosityMapping();
    extrapolateLevelsetOutside(m_fluidSdf);
    afterTransfer();
    extrapolateLevelsetInside(m_fluidSdf);
    m_stats.endStage(AFTER_TRANSFER);

    updateMaterials();
    applyBodyForces();
    m_stats.endStage(GRID_UPDATE);

    m_pressureMatrix = getPressureProjectionMatrix();
    m_pressureSolver.compute(m_pressureMatrix);
    m_stats.endStage(DECOMPOSITION);

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

void NBFlipSolver::advect()
{
    FlipSolver::advect();

    Vertex offsetU(0.f,0.5f);
    Vertex offsetV(0.5f,0.f);
    Vertex offsetCentered(0.5f,0.5f);

    Grid2d<float> advectedU(m_sizeI + 1, m_sizeJ, 0.f, OOBStrategy::OOB_EXTEND, 0.f, offsetU);
    Grid2d<float> advectedV(m_sizeI, m_sizeJ + 1, 0.f, OOBStrategy::OOB_EXTEND, 0.f, offsetV);
    Grid2d<float> advectedViscosity(m_sizeI, m_sizeJ, 0.f, OOBStrategy::OOB_EXTEND, 0.f, offsetCentered);

    std::vector<Range> rangesU = ThreadPool::i()->splitRange(advectedU.data().size());
    std::vector<Range> rangesV = ThreadPool::i()->splitRange(advectedV.data().size());
    std::vector<Range> rangesSdf = ThreadPool::i()->splitRange(m_advectedSdf.data().size());

    for(Range r : rangesU)
    {
        ThreadPool::i()->enqueue(&NBFlipSolver::eulerAdvectionThread,this,r,
                                 offsetU,std::ref(m_fluidVelocityGrid.velocityGridU()), std::ref(advectedU));
    }
    for(Range r : rangesV)
    {
        ThreadPool::i()->enqueue(&NBFlipSolver::eulerAdvectionThread,this,r,
                                 offsetV,std::ref(m_fluidVelocityGrid.velocityGridV()), std::ref(advectedV));
    }
    for(Range r : rangesSdf)
    {
        ThreadPool::i()->enqueue(&NBFlipSolver::eulerAdvectionThread,this,r,
                                 offsetCentered,std::ref(m_fluidSdf), std::ref(m_advectedSdf));
    }
    for(Range r : rangesSdf)
    {
        ThreadPool::i()->enqueue(&NBFlipSolver::eulerAdvectionThread,this,r,
                                 offsetCentered,std::ref(m_viscosityGrid), std::ref(advectedViscosity));
    }
    ThreadPool::i()->wait();

    m_viscosityGrid.swap(advectedViscosity);

    m_advectedVelocity.velocityGridU() = advectedU;
    m_advectedVelocity.velocityGridV() = advectedV;
}

void NBFlipSolver::afterTransfer()
{
    FlipSolver::afterTransfer();
    updateSdfFromSources();
    combineAdvectedGrids();
}

void NBFlipSolver::particleVelocityToGrid()
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
        ThreadPool::i()->enqueue(&NBFlipSolver::particleVelocityToGridThread,this,range,
                                 std::ref(uWeights), std::ref(vWeights));
    }
    ThreadPool::i()->wait();
}

void NBFlipSolver::particleVelocityToGridThread(Range r, Grid2d<float> &uWeights, Grid2d<float> &vWeights)
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
//                    if(m_fluidSdf.lerpolateAt(position) < m_combinationBand)
//                    {
//                        continue;
//                    }
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

void NBFlipSolver::reseedParticles()
{
    for (int pIndex = 0; pIndex < m_markerParticles.particleCount(); pIndex++)
    {
        Vertex pos = m_markerParticles.particlePosition(pIndex);
        int i = pos.x();
        int j = pos.y();
        float sdf = m_fluidSdf.lerpolateAt(pos);
        if(!m_materialGrid.isSource(i,j) && sdf < m_narrowBand)
        {
            m_markerParticles.markForDeath(pIndex);
            m_fluidParticleCounts.at(i,j) -= 1;
            continue;
        }

        if((m_materialGrid.isSource(i,j) || sdf < m_resamplingBand)
            && m_fluidParticleCounts.at(i,j) > 2*m_particlesPerCell)
        {
            m_markerParticles.markForDeath(pIndex);
            m_fluidParticleCounts.at(i,j) -= 1;
            continue;
        }
    }

    std::vector<float>& particleViscosities = m_markerParticles.particleProperties<float>
                                              (m_viscosityPropertyIndex);

    for (int i = 0; i < m_sizeI; i++)
    {
        for (int j = 0; j < m_sizeJ; j++)
        {
            if(m_materialGrid.isSource(i,j) || (m_fluidSdf.at(i,j) < m_resamplingBand
                                                  && m_fluidSdf.at(i,j) > m_narrowBand))
            {
                int particleCount = m_fluidParticleCounts.at(i,j);
                if(particleCount > 20)
                {
                    std::cout << "too many particles " << particleCount << " at " << i << ' ' << j;
                }
                int additionalParticles = m_particlesPerCell - particleCount;
                if(additionalParticles <= 0)
                {
                    continue;
                }
                for(int p = 0; p < additionalParticles; p++)
                {
                    Vertex pos = jitteredPosInCell(i,j);
                    float newSdf = m_fluidSdf.lerpolateAt(pos);
                    if((m_materialGrid.isSource(i,j) && m_fluidSdf.at(i,j) > m_resamplingBand)
                        || newSdf < m_narrowBand)
                    {
                        continue;
                    }
                    Vertex velocity = m_fluidVelocityGrid.velocityAt(pos);
                    //Vertex velocity = Vertex();
                    int emitterId = m_emitterId.at(i,j);
                    float viscosity = 0.f;
                    if(emitterId != -1)
                    {
                        viscosity = m_sources[emitterId].viscosity();
                    }
                    else
                    {
                        viscosity = m_viscosityGrid.interpolateAt(pos);
                    }
                    //float conc = m_sources[emitterId].concentrartion();
                    //float temp = m_sources[emitterId].temperature();
                    size_t pIdx = m_markerParticles.addMarkerParticle(pos,velocity);
                    particleViscosities[pIdx] = viscosity;
                }
            }
        }
    }
}

void NBFlipSolver::firstFrameInit()
{
    fluidSdfFromInitialFluid();
    m_fluidParticleCounts.fill(0);
    initialFluidSeed();
}

void NBFlipSolver::eulerAdvectionThread(Range range, Vertex offset, const Grid2d<float> &inputGrid, Grid2d<float> &outputGrid)
{
    std::vector<float>& dataOut = outputGrid.data();
    for(int idx = range.start; idx < range.end; idx++)
    {
        Index2d i2d = outputGrid.index2d(idx);
        Vertex currentPos = Vertex(i2d.i, i2d.j);
        Vertex prevPos = inverseRk4Integrate(currentPos,m_fluidVelocityGrid);
        dataOut[idx] = inputGrid.interpolateAt(prevPos);
    }
}

void NBFlipSolver::initialFluidSeed()
{
    std::vector<float>& particleViscosities = m_markerParticles.particleProperties<float>
                                              (m_viscosityPropertyIndex);

    for (int i = 0; i < m_sizeI; i++)
    {
        for (int j = 0; j < m_sizeJ; j++)
        {
            //if(m_fluidSdf.at(i,j) < 0.f)
            if(m_fluidSdf.at(i,j) < 0.f)
            {
                int particleCount = m_fluidParticleCounts.at(i,j);
                if(particleCount > 20)
                {
                    std::cout << "too many particles " << particleCount << " at " << i << ' ' << j;
                }
                int additionalParticles = m_particlesPerCell - particleCount;
                if(additionalParticles <= 0)
                {
                    continue;
                }
                for(int p = 0; p < additionalParticles; p++)
                {
                    Vertex pos = jitteredPosInCell(i,j);
                    float newSdf = m_fluidSdf.lerpolateAt(pos);
                    if(newSdf > m_resamplingBand || newSdf < m_narrowBand)
                    {
                        continue;
                    }
                    Vertex velocity = m_fluidVelocityGrid.velocityAt(pos);
                    //Vertex velocity = Vertex();
                    float viscosity = 0.f;
                    viscosity = m_viscosityGrid.interpolateAt(pos);
                    //float conc = m_sources[emitterId].concentrartion();
                    //float temp = m_sources[emitterId].temperature();
                    size_t pIdx = m_markerParticles.addMarkerParticle(pos,velocity);
                    particleViscosities[pIdx] = viscosity;
                }
            }
        }
    }
}

void NBFlipSolver::fluidSdfFromInitialFluid()
{
    for (int i = 0; i < m_sizeI; i++)
    {
        for (int j = 0; j < m_sizeJ; j++)
        {
            float dx = m_dx;
            float dist = std::numeric_limits<float>::max();
            int fluidId = -1;
            for(int fluidIdx = 0; fluidIdx < m_initialFluid.size(); fluidIdx++)
            {
                float sdf = m_initialFluid[fluidIdx].geometry().signedDistance(
                            (static_cast<float>(i)+0.5)*dx,
                            (static_cast<float>(j)+0.5)*dx) / dx;
                if(sdf < dist)
                {
                    dist = sdf;
                    fluidId = fluidIdx;
                }
            }
            m_fluidSdf.at(i,j) = dist;
            if(dist < 0)
            {
                m_materialGrid.at(i,j) = FluidMaterial::FLUID;
                if(fluidId != -1)
                {
                    m_viscosityGrid.at(i,j) = m_initialFluid[fluidId].viscosity();
                }
            }
        }
    }
}

void NBFlipSolver::updateSdfFromSources()
{
    for (int i = 0; i < m_sizeI; i++)
    {
        for (int j = 0; j < m_sizeJ; j++)
        {
            float dx = m_dx;
            float dist = std::numeric_limits<float>::max();
            int sourceId = -1;
            for(int sourceIdx = 0; sourceIdx < m_sources.size(); sourceIdx++)
            {
                float sdf = m_sources[sourceIdx].geometry().signedDistance(
                            (static_cast<float>(i))*dx,
                            (static_cast<float>(j))*dx) / dx;
                if(sdf < dist)
                {
                    dist = sdf;
                    sourceId = sourceIdx;
                }
            }
            m_fluidSdf.at(i,j) = std::min(dist,m_fluidSdf.at(i,j));
            if(dist < 0)
            {
                //m_materialGrid.at(i,j) = FluidMaterial::SOURCE;
                if(sourceId != -1)
                {
                    m_viscosityGrid.at(i,j) = m_sources[sourceId].viscosity();
                }
            }
        }
    }
}

void NBFlipSolver::combineAdvectedGrids()
{
    combineLevelset();
    combineVelocityGrid();
}

void NBFlipSolver::combineLevelset()
{
    const float h = 1.f;
    for (int i = 0; i < m_sizeI; i++)
    {
        for (int j = 0; j < m_sizeJ; j++)
        {
            //m_fluidSdf.at(i,j) = m_advectedSdf.at(i,j);
            m_fluidSdf.at(i,j) = std::min(m_advectedSdf.at(i,j) + h, m_fluidSdf.at(i,j));
        }
    }
}

void NBFlipSolver::combineVelocityGrid()
{
    auto nbcombine = [this](float vAdvected, float vParticle, float sdf)
    {
        if(sdf > m_combinationBand)
        {
            return vParticle;
        }
        else
        {
            return vAdvected;
        }
    };
    for (int i = 0; i < m_sizeI + 1; i++)
    {
        for (int j = 0; j < m_sizeJ; j++)
        {
            Vertex samplePosition = Vertex(i, 0.5f + j);
            float advectedVelocity = m_advectedVelocity.getU(i,j);
            float particleVelocity = m_fluidVelocityGrid.getU(i,j);
            float sdf = m_fluidSdf.lerpolateAt(samplePosition);
            m_fluidVelocityGrid.setU(i,j,nbcombine(advectedVelocity, particleVelocity, sdf));
        }
    }

    for (int i = 0; i < m_sizeI; i++)
    {
        for (int j = 0; j < m_sizeJ + 1; j++)
        {
            Vertex samplePosition = Vertex(0.5f + i, j);
            float advectedVelocity = m_advectedVelocity.getV(i,j);
            float particleVelocity = m_fluidVelocityGrid.getV(i,j);
            float sdf = m_fluidSdf.lerpolateAt(samplePosition);
            m_fluidVelocityGrid.setV(i,j,nbcombine(advectedVelocity, particleVelocity, sdf));
        }
    }
}

Vertex NBFlipSolver::inverseRk4Integrate(Vertex newPosition, StaggeredVelocityGrid &grid)
{
    float factor = -m_stepDt;
    Vertex k1 = factor*grid.velocityAt(newPosition);
    Vertex k2 = factor*grid.velocityAt(newPosition + 0.5f*k1);
    Vertex k3 = factor*grid.velocityAt(newPosition + 0.5f*k2);
    Vertex k4 = factor*grid.velocityAt(newPosition + k3);

    return newPosition + (1.0f/6.0f)*(k1 + 2.f*k2 + 2.f*k3 + k4);
}

const StaggeredVelocityGrid &NBFlipSolver::advectedVelocityGrid() const
{
    return m_advectedVelocity;
}

void NBFlipSolver::centeredParamsToGrid()
{
    Grid2d<float> centeredWeights(m_sizeI,m_sizeJ,1e-10f);

    //m_viscosityGrid.fill(0.f);
    m_divergenceControl.fill(0.f);
    m_knownCenteredParams.fill(false);

    std::vector<Range> ranges = ThreadPool::i()->splitRange(m_viscosityGrid.linearSize(),1,8);

    for(Range& range : ranges)
    {
        ThreadPool::i()->enqueue(&NBFlipSolver::centeredParamsToGridThread,this,range,
                                 std::ref(centeredWeights));
    }
    ThreadPool::i()->wait();
}

void NBFlipSolver::centeredParamsToGridThread(Range r, Grid2d<float> &cWeights)
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
                        if(!m_knownCenteredParams.at(i2d.i,i2d.j))
                        {
                            m_knownCenteredParams.at(i2d.i,i2d.j) = true;
                            m_viscosityGrid.at(i2d) = 0.f;
                        }
                        m_viscosityGrid.at(i2d.i,i2d.j) += weightCentered * particleViscosities.at(particleIdx);
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
