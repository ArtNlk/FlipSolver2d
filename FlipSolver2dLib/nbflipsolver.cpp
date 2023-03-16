#include "nbflipsolver.h"
#include "flipsolver2d.h"
#include "sdfgrid.h"
#include "simsettings.h"
#include "staggeredvelocitygrid.h"
#include <algorithm>

NBFlipSolver::NBFlipSolver(int sizeI, int sizeJ, int extrapRadius, bool vonNeumannNeighbors):
    FlipSolver(sizeI, sizeJ, extrapRadius, vonNeumannNeighbors),
    m_advectedVelocity(sizeI, sizeJ),
    m_advectedSdf(sizeI, sizeJ),
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

    Grid2d<float> advectedU(m_sizeI + 1, m_sizeJ, 0.f, OOBStrategy::OOB_EXTEND, 0.f, Vertex(0.5f, 0.f));
    Grid2d<float> advectedV(m_sizeI, m_sizeJ + 1, 0.f, OOBStrategy::OOB_EXTEND, 0.f, Vertex(0.f, 0.5f));

    for(int i = 0; i < m_sizeI + 1; i++)
    {
        for(int j = 0; j < m_sizeJ; j++)
        {
            Vertex currentPos(0.5f + i,j);
            Vertex prevPos = inverseRk3Integrate(currentPos,m_fluidVelocityGrid);
            advectedU.at(i,j) = m_fluidVelocityGrid.velocityGridU().interpolateAt(prevPos);
        }
    }

    for(int i = 0; i < m_sizeI; i++)
    {
        for(int j = 0; j < m_sizeJ + 1; j++)
        {
            Vertex currentPos(i,0.5f + j);
            Vertex prevPos = inverseRk3Integrate(currentPos,m_fluidVelocityGrid);
            advectedV.at(i,j) = m_fluidVelocityGrid.velocityGridV().interpolateAt(prevPos);
        }
    }

    for(int i = 0; i < m_sizeI; i++)
    {
        for(int j = 0; j < m_sizeJ; j++)
        {
            Vertex currentPos(0.5f + i,0.5f + j);
            Vertex prevPos = inverseRk3Integrate(currentPos,m_fluidVelocityGrid);
            m_advectedSdf.at(i,j) = m_fluidSdf.interpolateAt(prevPos);
        }
    }
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
        if(sdf < -3.f)
        {
            m_markerParticles.erase(m_markerParticles.cbegin() + pIndex);
            m_fluidParticleCounts.at(i,j) -= 1;
            pIndex--;
            continue;
        }

        if(sdf < -1.f && m_fluidParticleCounts.at(i,j) > 2*SimSettings::particlesPerCell())
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
            if((m_materialGrid.isSource(i,j) || m_fluidSdf.at(i,j) < -1.f) && m_fluidSdf.at(i,j) > -3.f)
            {
                int particleCount = m_fluidParticleCounts.at(i,j);
                if(particleCount > 20)
                {
                    std::cout << "too many particles " << particleCount << " at " << i << ' ' << j;
                }
                int additionalParticles = SimSettings::particlesPerCell() - particleCount;
                if(additionalParticles <= 0)
                {
                    continue;
                }
                for(int p = 0; p < additionalParticles; p++)
                {
                    Vertex pos = jitteredPosInCell(i,j);
                    float newSdf = m_fluidSdf.interpolateAt(pos);
                    if((!m_materialGrid.isSource(i,j) && m_fluidSdf.at(i,j) > -1.f)
                        || newSdf < -3.f)
                    {
                        continue;
                    }
                    Vertex velocity = m_fluidVelocityGrid.velocityAt(pos);
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
                int additionalParticles = SimSettings::particlesPerCell() - particleCount;
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
            float dx = SimSettings::dx();
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
            float dx = SimSettings::dx();
            float dist = std::numeric_limits<float>::max();
            int sourceId = -1;
            for(int sourceIdx = 0; sourceIdx < m_sources.size(); sourceIdx++)
            {
                float sdf = m_sources[sourceIdx].geometry().signedDistance(
                            (static_cast<float>(i)+0.5)*dx,
                            (static_cast<float>(j)+0.5)*dx) / dx;
                if(sdf < dist)
                {
                    dist = sdf;
                    sourceId = sourceIdx;
                }
            }
            m_fluidSdf.at(i,j) = std::min(dist,m_fluidSdf.at(i,j));
            if(dist < 0)
            {
                m_materialGrid.at(i,j) = FluidMaterial::FLUID;
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

Vertex NBFlipSolver::inverseRk3Integrate(Vertex newPosition, StaggeredVelocityGrid &grid)
{
    float dt = SimSettings::stepDt();
    Vertex k1 = grid.velocityAt(newPosition);
    Vertex k2 = grid.velocityAt(newPosition - 1.f/2.f * dt * k1);
    Vertex k3 = grid.velocityAt(newPosition - 3.f/4.f * dt * k2);

    return newPosition - (2.f/9.f) * dt * k1 - (3.f/9.f) * dt * k2 - (4.f/9.f) * dt * k3;
}
