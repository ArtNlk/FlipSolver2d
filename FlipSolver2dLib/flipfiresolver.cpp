#include "flipfiresolver.h"
#include "flipsmokesolver.h"
#include "grid2d.h"


FlipFireSolver::FlipFireSolver(const FireSolverParameters *p):
    FlipSmokeSolver(p),
    m_fuel(p->gridSizeI, p->gridSizeJ, 0.f, OOBStrategy::OOB_CONST, 0.f),
    m_oxidizer(p->gridSizeI, p->gridSizeJ, 0.f, OOBStrategy::OOB_CONST, 0.f),
    m_ignitionTemperature(p->ignitionTemperature),
    m_burnRate(p->burnRate),
    m_smokeProportion(p->smokeProportion),
    m_heatProportion(p->heatProportion),
    m_divergenceProportion(p->divergenceProportion)
{
}

void FlipFireSolver::afterTransfer()
{
    FlipSmokeSolver::afterTransfer();
    for (int i = 0; i < m_sizeI; i++)
    {
        for (int j = 0; j < m_sizeJ; j++)
        {
            if(m_materialGrid.isSource(i,j))
            {
                int emitterId = m_emitterId.at(i,j);
                m_fuel.at(i,j) = m_sources[emitterId].fuel();
            }
        }
    }
//    combustionUpdate();
}

void FlipFireSolver::combustionUpdate()
{
    std::vector<Range> ranges = ThreadPool::i()->splitRange(m_markerParticles.particleCount());

    for(Range& range : ranges)
    {
        ThreadPool::i()->enqueue(&FlipFireSolver::combustionUpdateThread,this,range);
    }
    ThreadPool::i()->wait();
}

void FlipFireSolver::combustionUpdateThread(Range range)
{
    std::vector<float>& particleTemperatures = m_markerParticles.particleProperties<float>
                                               (m_temperatureIndex);
    std::vector<float>& particleConcentrations = m_markerParticles.particleProperties<float>
                                                 (m_concentrationIndex);
    std::vector<float>& particleFuels = m_markerParticles.particleProperties<float>
                                        (m_fuelPropertyIndex);
    std::vector<float>& testValues = m_markerParticles.particleProperties<float>
                                        (m_testValuePropertyIndex);

    for(int i = range.start; i < range.end; i++)
    {
        if(particleTemperatures.at(i) > m_ignitionTemperature && particleFuels.at(i) > 0.f)
        {
            float burntFuel = std::min(m_stepDt * m_burnRate,particleFuels.at(i));
            particleFuels.at(i) -= burntFuel;
            testValues.at(i) = particleFuels.at(i);
            particleConcentrations.at(i) += m_smokeProportion * burntFuel;
            particleTemperatures.at(i) += m_heatProportion * burntFuel;
            //m_divergenceControl.data().at(i) -= m_divergenceProportion * burntFuel;
        }
    }
}

void FlipFireSolver::centeredParamsToGrid()
{
    Grid2d<float> centeredWeights(m_sizeI,m_sizeJ,1e-10f);
    m_fuel.fill(0.f);
    m_smokeConcentration.fill(0.f);
    m_temperature.fill(0.f);
    m_divergenceControl.fill(0.f);
    m_knownCenteredParams.fill(false);

    std::vector<float>& particleTemperatures = m_markerParticles.particleProperties<float>
                                               (m_temperatureIndex);

    std::vector<float>& particleConcentrations = m_markerParticles.particleProperties<float>
                                                 (m_concentrationIndex);

    std::vector<float>& particleFuels = m_markerParticles.particleProperties<float>
                                                 (m_fuelPropertyIndex);

    for(size_t idx = 0; idx < m_fluidVelocityGrid.linearSize(); idx++)
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
                        centeredWeights.at(i2d.i,i2d.j) += weightCentered;
                        m_temperature.at(i2d.i,i2d.j) +=
                            weightCentered * particleTemperatures.at(particleIdx);

                        m_smokeConcentration.at(i2d.i,i2d.j) +=
                            weightCentered * particleConcentrations.at(particleIdx);

                        m_fuel.at(i2d.i,i2d.j) +=
                            weightCentered * particleFuels.at(particleIdx);

                        m_knownCenteredParams.at(i2d.i,i2d.j) = true;
                    }
                }
            }
        }
    }

    for (int i = 0; i < m_sizeI + 1; i++)
    {
        for (int j = 0; j < m_sizeJ + 1; j++)
        {
            if(centeredWeights.inBounds(i,j))
            {
                if(m_knownCenteredParams.at(i,j))
                {
                    m_temperature.at(i,j) /= centeredWeights.at(i,j);
                    m_smokeConcentration.at(i,j) /= centeredWeights.at(i,j);
                    m_fuel.at(i,j) /= centeredWeights.at(i,j);
                }
                else
                {
                    m_temperature.at(i,j) = m_ambientTemperature;
                }
                //m_testGrid.at(i,j) = m_temperature.at(i,j) / 1000.f;
            }
        }
    }
}

void FlipFireSolver::reseedParticles()
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

    std::vector<float>& particleTemperatures = m_markerParticles.particleProperties<float>
                                               (m_temperatureIndex);

    std::vector<float>& particleConcentrations = m_markerParticles.particleProperties<float>
                                                 (m_concentrationIndex);

    std::vector<float>& particleFuels = m_markerParticles.particleProperties<float>
                                        (m_fuelPropertyIndex);

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
                    float conc = m_smokeConcentration.interpolateAt(pos);
                    float temp = m_temperature.interpolateAt(pos);
                    float fuel = m_fuel.interpolateAt(pos);
                    size_t pIdx = m_markerParticles.addMarkerParticle(pos,velocity);
                    particleTemperatures[pIdx] = temp;
                    particleConcentrations[pIdx] = conc;
                    particleFuels[pIdx] = fuel;
                    //testValues[pIdx] = temp / 1000.f;
                }
            }
        }
    }
}

void FlipFireSolver::particleUpdate()
{
    FlipSmokeSolver::particleUpdate();
    combustionUpdate();
//    for(int i = m_markerParticles.size() - 1; i >= 0; i--)
//    {
//        MarkerParticle &p = m_markerParticles[i];
//        p.fuel = m_fuel.lerpolateAt(p.position);
    //    }
}

void FlipFireSolver::initAdditionalParameters()
{
    FlipSmokeSolver::initAdditionalParameters();

    m_fuelPropertyIndex = m_markerParticles.addParticleProperty<float>();
}
