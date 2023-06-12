#include "nbflipsolver.h"
#include "flipsolver2d.h"
#include "sdfgrid.h"

#include "staggeredvelocitygrid.h"
#include <algorithm>

NBFlipSolver::NBFlipSolver(const NBFlipParameters *p):
    FlipSolver(p),
    m_advectedVelocity(p->gridSizeI, p->gridSizeJ),
    m_advectedSdf(p->gridSizeI, p->gridSizeJ),
    m_uWeights(m_fluidVelocityGrid.velocityGridU().sizeI(),
            m_fluidVelocityGrid.velocityGridU().sizeJ(),1e-10f),
    m_vWeights(m_fluidVelocityGrid.velocityGridV().sizeI(),
            m_fluidVelocityGrid.velocityGridV().sizeJ(),1e-10f)
{
}

void NBFlipSolver::step()
{
    advect();
    particleToGrid();
    m_fluidVelocityGrid.extrapolate(10);
    m_savedFluidVelocityGrid = m_fluidVelocityGrid;
    updateSdf();
    updateLinearFluidViscosityMapping();
    afterTransfer();
    extrapolateLevelsetInside(m_fluidSdf);
    updateMaterials();
    applyBodyForces();
    project();
    updateVelocityFromSolids();
    applyViscosity();
    project();
    //m_fluidVelocityGrid.extrapolate(10);
    particleUpdate();
    countParticles();
    reseedParticles();
}

void NBFlipSolver::advect()
{
    FlipSolver::advect();

    Vertex offsetU(0.5f,0.f);
    Vertex offsetV(0.f,0.5f);
    Vertex offsetCentered(0.5f,0.5f);

    Grid2d<float> advectedU(m_sizeI + 1, m_sizeJ, 0.f, OOBStrategy::OOB_EXTEND, 0.f, offsetU);
    Grid2d<float> advectedV(m_sizeI, m_sizeJ + 1, 0.f, OOBStrategy::OOB_EXTEND, 0.f, offsetV);

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
    ThreadPool::i()->wait();

    m_advectedVelocity.velocityGridU() = advectedU;
    m_advectedVelocity.velocityGridV() = advectedV;
}

void NBFlipSolver::afterTransfer()
{
    //FlipSolver::afterTransfer();
    updateSdfFromSources();
    combineLevelset();
    combineVelocityGrid();
}

void NBFlipSolver::particleVelocityToGrid()
{
    m_uWeights.fill(1e-10f);
    m_vWeights.fill(1e-10f);
    m_fluidVelocityGrid.velocityGridU().fill(0.f);
    m_fluidVelocityGrid.velocityGridV().fill(0.f);
    m_fluidVelocityGrid.uSampleValidityGrid().fill(false);
    m_fluidVelocityGrid.vSampleValidityGrid().fill(false);

    for(MarkerParticle &p : m_markerParticles)
    {
        int i = simmath::integr(p.position.x());
        int j = simmath::integr(p.position.y());
        //Run over all cells that this particle might affect
        for (int iOffset = -3; iOffset < 3; iOffset++)
        {
            for (int jOffset = -3; jOffset < 3; jOffset++)
            {
                int iIdx = i+iOffset;
                int jIdx = j+jOffset;
                float weightU = simmath::quadraticBSpline(p.position.x() - (iIdx),
                                                     p.position.y() - (static_cast<float>(jIdx) + 0.5f));

                float weightV = simmath::quadraticBSpline(p.position.x() - (static_cast<float>(iIdx) + 0.5f),
                                                     p.position.y() - (jIdx));
                if(m_uWeights.inBounds(iIdx,jIdx))
                {
                    if(std::abs(weightU) > 1e-9f && std::abs(weightV) > 1e-9f)
                    {
                        m_uWeights.at(iIdx,jIdx) += weightU;
                        m_fluidVelocityGrid.u(iIdx,jIdx) += weightU * (p.velocity.x());
                        m_fluidVelocityGrid.setUValidity(iIdx,jIdx,true);
                    }
                }

                if(m_vWeights.inBounds(iIdx,jIdx))
                {
                    if(std::abs(weightU) > 1e-9f && std::abs(weightV) > 1e-9f)
                    {
                        m_vWeights.at(iIdx,jIdx) += weightV;
                        m_fluidVelocityGrid.v(iIdx,jIdx) += weightV * (p.velocity.y());
                        m_fluidVelocityGrid.setVValidity(iIdx,jIdx,true);
                    }
                }
            }
        }
    }

    for (int i = 0; i < m_sizeI + 1; i++)
    {
        for (int j = 0; j < m_sizeJ + 1; j++)
        {
            if(m_fluidVelocityGrid.velocityGridU().inBounds(i,j))
            {
                m_fluidVelocityGrid.u(i,j) /= m_uWeights.at(i,j);
            }
            if(m_fluidVelocityGrid.velocityGridV().inBounds(i,j))
            {
                m_fluidVelocityGrid.v(i,j) /= m_vWeights.at(i,j);
            }
        }
    }
}

void NBFlipSolver::reseedParticles()
{
    for (int pIndex = 0; pIndex < m_markerParticles.size(); pIndex++)
    {
        float sdf = m_fluidSdf.interpolateAt(m_markerParticles[pIndex].position);
        int i = m_markerParticles[pIndex].position.x();
        int j = m_markerParticles[pIndex].position.y();
        if(!m_materialGrid.isSource(i,j) && sdf < -3.f)
        {
            m_markerParticles.erase(m_markerParticles.cbegin() + pIndex);
            m_fluidParticleCounts.at(i,j) -= 1;
            pIndex--;
            continue;
        }

        if((m_materialGrid.isSource(i,j) || sdf < -1.f) && m_fluidParticleCounts.at(i,j) > 2*m_particlesPerCell)
        {
            m_markerParticles.erase(m_markerParticles.cbegin() + pIndex);
            m_fluidParticleCounts.at(i,j) -= 1;
            pIndex--;
            continue;
        }
    }

    for (int i = 0; i < m_sizeI; i++)
    {
        for (int j = 0; j < m_sizeJ; j++)
        {
            //if(m_fluidSdf.at(i,j) < 0.f)
            if(m_materialGrid.isSource(i,j) || (m_fluidSdf.at(i,j) < -1.f && m_fluidSdf.at(i,j) > -3.f))
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
                    float newSdf = m_fluidSdf.interpolateAt(pos);
                    if((m_materialGrid.isSource(i,j) && m_fluidSdf.at(i,j) > -1.f)
                        || newSdf < -3.f)
                    {
                        continue;
                    }
                    Vertex velocity = m_materialGrid.isSource(i,j) ?
                                          Vertex() : m_fluidVelocityGrid.velocityAt(pos);
                    float viscosity = m_viscosityGrid.interpolateAt(pos);
                    addMarkerParticle(MarkerParticle{pos,velocity,viscosity});
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
        Vertex currentPos = Vertex(i2d.m_i, i2d.m_j) + offset;
        Vertex prevPos = inverseRk4Integrate(currentPos,m_fluidVelocityGrid);
        dataOut[idx] = inputGrid.interpolateAt(prevPos);
    }
}

void NBFlipSolver::initialFluidSeed()
{
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
//                    float newSdf = m_fluidSdf.interpolateAt(pos);
//                    if(newSdf > -1.f || newSdf < -3.f)
//                    {
//                        continue;
//                    }
                    Vertex velocity = m_fluidVelocityGrid.velocityAt(pos);
                    float viscosity = m_viscosityGrid.interpolateAt(pos);
                    addMarkerParticle(MarkerParticle{pos,velocity,viscosity});
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
    float h = 1.f;
    for (int i = 0; i < m_sizeI; i++)
    {
        for (int j = 0; j < m_sizeJ; j++)
        {
            m_fluidSdf.at(i,j) = std::min(m_advectedSdf.at(i,j) + h, m_fluidSdf.at(i,j));
        }
    }
}

void NBFlipSolver::combineVelocityGrid()
{
    auto nbcombine = [](float vAdvected, float vParticle, float sdf)
    {
        const float r = 2.f;
        if(sdf >= -r)
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
            Vertex samplePosition = Vertex(0.5f + i, j);
            float advectedVelocity = m_advectedVelocity.getU(i,j);
            float particleVelocity = m_fluidVelocityGrid.getU(i,j);
            float sdf = m_fluidSdf.interpolateAt(samplePosition);
            m_fluidVelocityGrid.setU(i,j,nbcombine(advectedVelocity, particleVelocity, sdf));
        }
    }

    for (int i = 0; i < m_sizeI; i++)
    {
        for (int j = 0; j < m_sizeJ + 1; j++)
        {
            Vertex samplePosition = Vertex(i, 0.5f + j);
            float advectedVelocity = m_advectedVelocity.getV(i,j);
            float particleVelocity = m_fluidVelocityGrid.getV(i,j);
            float sdf = m_fluidSdf.interpolateAt(samplePosition);
            m_fluidVelocityGrid.setV(i,j,nbcombine(advectedVelocity, particleVelocity, sdf));
        }
    }
}

Vertex NBFlipSolver::inverseRk4Integrate(Vertex newPosition, StaggeredVelocityGrid &grid)
{
    float factor = m_stepDt / m_dx;
    Vertex k1 = factor*grid.velocityAt(newPosition);
    Vertex k2 = factor*grid.velocityAt(newPosition - 0.5f*k1);
    Vertex k3 = factor*grid.velocityAt(newPosition - 0.5f*k2);
    Vertex k4 = factor*grid.velocityAt(newPosition - k3);

    return newPosition - (1.0f/6.0f)*(k1 - 2.f*k2 - 2.f*k3 - k4);
}
