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
    combustionUpdate();
}

void FlipFireSolver::combustionUpdate()
{
    std::vector<Range> ranges = ThreadPool::i()->splitRange(m_fuel.data().size());

    for(Range& range : ranges)
    {
        ThreadPool::i()->enqueue(&FlipFireSolver::combustionUpdateThread,this,range);
    }
    ThreadPool::i()->wait();
}

void FlipFireSolver::combustionUpdateThread(Range range)
{
    for(int i = range.start; i < range.end; i++)
    {
        if(m_temperature.data().at(i) > m_ignitionTemperature && m_fuel.data().at(i) > 0.f)
        {
            float burntFuel = std::min(m_stepDt * m_burnRate,m_fuel.data().at(i));
            m_fuel.data().at(i) -= burntFuel;
            m_smokeConcentration.data().at(i) += m_smokeProportion * burntFuel;
            m_temperature.data().at(i) += m_heatProportion * burntFuel;
            m_divergenceControl.data().at(i) -= m_divergenceProportion * burntFuel;
        }
    }
}

void FlipFireSolver::centeredParamsToGrid()
{
    FlipSmokeSolver::centeredParamsToGrid();
    Grid2d<float> centeredWeights(m_sizeI,m_sizeJ,1e-10f);

    m_fuel.fill(0.f);
    m_knownCenteredParams.fill(false);

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
                if(inBounds(iIdx,jIdx) && !m_materialGrid.isSource(i,j))
                {
                    float weightCentered = simmath::quadraticBSpline(p.position.x() - (iIdx + 0.5f),
                                                         p.position.y() - (jIdx + 0.5f));
                    if(std::abs(weightCentered) > 1e-6f)
                    {
                        centeredWeights.at(iIdx,jIdx) += weightCentered;
                        m_fuel.at(iIdx,jIdx) += weightCentered * p.fuel;
                        m_knownCenteredParams.at(iIdx,jIdx) = true;
                    }
                }
            }
        }
    }

    for (int i = 0; i < m_sizeI + 1; i++)
    {
        for (int j = 0; j < m_sizeJ + 1; j++)
        {
            if(centeredWeights.inBounds(i,j) && !m_materialGrid.isSource(i,j))
            {
                if(m_knownCenteredParams.at(i,j))
                {
                    m_fuel.at(i,j) /= centeredWeights.at(i,j);
                }
            }
        }
    }
}

void FlipFireSolver::reseedParticles()
{
    for (int pIndex = 0; pIndex < m_markerParticles.size(); pIndex++)
    {
        int i = m_markerParticles[pIndex].position.x();
        int j = m_markerParticles[pIndex].position.y();

        if(m_fluidParticleCounts.at(i,j) > 2*m_particlesPerCell)
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
            int particleCount = m_fluidParticleCounts.at(i,j);
            if(particleCount > 20)
            {
                std::cout << "too many particles " << particleCount << " at " << i << ' ' << j;
            }
            //std::cout << particleCount << " at " << i << " , " << j << std::endl;
            int additionalParticles = m_particlesPerCell - particleCount;
            if(additionalParticles <= 0)
            {
                continue;
            }
            if(m_materialGrid.isSource(i,j))
            {
                for(int p = 0; p < additionalParticles; p++)
                {
                    Vertex pos = jitteredPosInCell(i,j);
                    //Vertex velocity = m_fluidVelocityGrid.velocityAt(pos);
                    Vertex velocity = Vertex();
                    int emitterId = m_emitterId.at(i,j);
                    float viscosity = m_sources[emitterId].viscosity();
                    float conc = m_sources[emitterId].concentrartion();
                    float temp = m_sources[emitterId].temperature();
                    float fuel = m_sources[emitterId].fuel();
                    addMarkerParticle(MarkerParticle{pos,velocity,viscosity,temp,conc,fuel});
                }
            }
        }
    }
}

void FlipFireSolver::particleUpdate()
{
    FlipSmokeSolver::particleUpdate();
    for(int i = m_markerParticles.size() - 1; i >= 0; i--)
    {
        MarkerParticle &p = m_markerParticles[i];
        p.fuel = m_fuel.lerpolateAt(p.position);
    }
}
