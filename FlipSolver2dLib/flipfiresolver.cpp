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
    for (ssize_t i = 0; i < m_sizeI; i++)
    {
        for (ssize_t j = 0; j < m_sizeJ; j++)
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
    for(size_t binIdx = 0; binIdx < m_markerParticles.bins().linearSize(); binIdx++)
    {
        ParticleBin& currentBin = m_markerParticles.bins().data().at(binIdx);
        std::vector<float>& particleTemperatures = currentBin.particleProperties<float>
                                                   (m_temperatureIndex);

        std::vector<float>& particleConcentrations = currentBin.particleProperties<float>
                                                     (m_concentrationIndex);

        std::vector<float>& particleFuels = currentBin.particleProperties<float>
                                                     (m_fuelPropertyIndex);

        for(size_t particleIdx = 0; particleIdx < currentBin.size(); particleIdx++)
        {
            if(particleTemperatures.at(particleIdx) > m_ignitionTemperature && particleFuels.at(particleIdx) > 0.f)
            {
                float burntFuel = std::min(m_stepDt * m_burnRate,particleFuels.at(particleIdx));
                particleFuels.at(particleIdx) -= burntFuel;
                //testValues.at(i) = particleFuels.at(i);
                particleConcentrations.at(particleIdx) += m_smokeProportion * burntFuel;
                particleTemperatures.at(particleIdx) += m_heatProportion * burntFuel;
                //m_divergenceControl.data().at(i) -= m_divergenceProportion * burntFuel;
            }
        }
    }
}

void FlipFireSolver::gridCombustionUpdateThread(Range range)
{
    std::vector<float>& temperatureData = m_temperature.data();
    std::vector<float>& concentrationData = m_smokeConcentration.data();
    std::vector<float>& fuelData = m_fuel.data();
    std::vector<float>& testGridData = m_testGrid.data();

    for(size_t i = range.start; i < range.end; i++)
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
    for (ssize_t i = 0; i < m_sizeI; i++)
    {
        for (ssize_t j = 0; j < m_sizeJ; j++)
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

            ParticleBin& bin = m_markerParticles.binForGridIdx(i,j);
            auto& tempVec = bin.particleProperties<float>(m_temperatureIndex);
            auto& concVec = bin.particleProperties<float>(m_concentrationIndex);
            auto& fuelVec = bin.particleProperties<float>(m_fuelPropertyIndex);

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
                    float conc = m_smokeConcentration.interpolateAt(pos);
                    float temp = m_temperature.interpolateAt(pos);
                    float fuel = m_fuel.interpolateAt(pos);
                    size_t pIdx = bin.addMarkerParticle(pos,velocity);
                    tempVec[pIdx] = temp;
                    concVec[pIdx] = conc;
                    fuelVec[pIdx] = fuel;
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
}

void FlipFireSolver::eulerAdvectParameters()
{
    FlipSmokeSolver::eulerAdvectParameters();

    Vec3 offsetCentered(0.5f,0.5f);
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
