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
    if(m_parameterHandlingMethod == PARTICLE)
    {
        std::vector<Range> ranges = ThreadPool::i()->splitRange(m_markerParticles.particleCount());

        for(Range& range : ranges)
        {
            ThreadPool::i()->enqueue(&FlipFireSolver::particleCombustionUpdateThread,this,range);
        }
    }
    else
    {
        std::vector<Range> ranges = ThreadPool::i()->splitRange(m_fuel.linearSize());

        for(Range& range : ranges)
        {
            ThreadPool::i()->enqueue(&FlipFireSolver::gridCombustionUpdateThread,this,range);
        }
    }
    ThreadPool::i()->wait();
}

void FlipFireSolver::particleCombustionUpdateThread(Range range)
{
    // std::vector<float>& particleTemperatures = m_markerParticles.particleProperties<float>
    //                                            (m_temperatureIndex);
    // std::vector<float>& particleConcentrations = m_markerParticles.particleProperties<float>
    //                                              (m_concentrationIndex);
    // std::vector<float>& particleFuels = m_markerParticles.particleProperties<float>
    //                                     (m_fuelPropertyIndex);
    // std::vector<float>& testValues = m_markerParticles.particleProperties<float>
    //                                     (m_testValuePropertyIndex);

    // for(int i = range.start; i < range.end; i++)
    // {
    //     if(particleTemperatures.at(i) > m_ignitionTemperature && particleFuels.at(i) > 0.f)
    //     {
    //         float burntFuel = std::min(m_stepDt * m_burnRate,particleFuels.at(i));
    //         particleFuels.at(i) -= burntFuel;
    //         testValues.at(i) = particleFuels.at(i);
    //         particleConcentrations.at(i) += m_smokeProportion * burntFuel;
    //         particleTemperatures.at(i) += m_heatProportion * burntFuel;
    //         //m_divergenceControl.data().at(i) -= m_divergenceProportion * burntFuel;
    //     }
    // }
}

void FlipFireSolver::gridCombustionUpdateThread(Range range)
{
    std::vector<float>& temperatureData = m_temperature.data();
    std::vector<float>& concentrationData = m_smokeConcentration.data();
    std::vector<float>& fuelData = m_fuel.data();
    std::vector<float>& testGridData = m_testGrid.data();

    for(int i = range.start; i < range.end; i++)
    {
        if(temperatureData.at(i) > m_ignitionTemperature && fuelData.at(i) > 0.f)
        {
            float burntFuel = std::min(m_stepDt * m_burnRate,fuelData.at(i));
            fuelData.at(i) -= burntFuel;
            testGridData.at(i) = fuelData.at(i);
            concentrationData.at(i) += m_smokeProportion * burntFuel;
            temperatureData.at(i) += m_heatProportion * burntFuel;
            //m_divergenceControl.data().at(i) -= m_divergenceProportion * burntFuel;
        }
    }
}

void FlipFireSolver::reseedParticles()
{
    // for (int pIndex = 0; pIndex < m_markerParticles.particleCount(); pIndex++)
    // {
    //     Vertex pos = m_markerParticles.particlePosition(pIndex);
    //     int i = pos.x();
    //     int j = pos.y();

    //     if(m_fluidParticleCounts.at(i,j) > 4*m_particlesPerCell)
    //     {
    //         // m_markerParticles.markForDeath(pIndex);
    //         m_fluidParticleCounts.at(i,j) -= 1;
    //         //pIndex--;
    //         continue;
    //     }
    // }

    // std::vector<float>& particleTemperatures = m_markerParticles.particleProperties<float>
    //                                            (m_temperatureIndex);

    // std::vector<float>& particleConcentrations = m_markerParticles.particleProperties<float>
    //                                              (m_concentrationIndex);

    // std::vector<float>& particleFuels = m_markerParticles.particleProperties<float>
    //                                     (m_fuelPropertyIndex);

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
                    // size_t pIdx = m_markerParticles.addMarkerParticle(pos,velocity);
                    // particleTemperatures[pIdx] = temp;
                    // particleConcentrations[pIdx] = conc;
                    // particleFuels[pIdx] = fuel;
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

void FlipFireSolver::eulerAdvectParameters()
{
    FlipSmokeSolver::eulerAdvectParameters();

    Vertex offsetCentered(0.5f,0.5f);
    Grid2d<float> advectedFuel(m_sizeI, m_sizeJ, 0.f, OOBStrategy::OOB_EXTEND, 0.f, offsetCentered);

    std::vector<Range> ranges = ThreadPool::i()->splitRange(advectedFuel.data().size());
    for(Range r : ranges)
    {
        ThreadPool::i()->enqueue(&FlipSolver::eulerAdvectionThread,this,r,
                                 std::ref(m_fuel),std::ref(advectedFuel));
    }
    ThreadPool::i()->wait();

    m_fuel = advectedFuel;
}

void FlipFireSolver::initAdditionalParameters()
{
    FlipSmokeSolver::initAdditionalParameters();

    m_fuelPropertyIndex = m_markerParticles.addParticleProperty<float>();
}
