#include "flipsmokesolver.h"

#include "flipsolver2d.h"
#include "mathfuncs.h"

#include <cmath>

FlipSmokeSolver::FlipSmokeSolver(const SmokeSolverParameters* p):
    FlipSolver(p),
    m_temperature(p->gridSizeI, p->gridSizeJ, p->ambientTemperature, OOBStrategy::OOB_CONST, p->ambientTemperature),
    m_smokeConcentration(p->gridSizeI, p->gridSizeJ, 0.f, OOBStrategy::OOB_CONST, 0.f),
    m_ambientTemperature(p->ambientTemperature),
    m_temperatureDecayRate(p->temperatureDecayRate),
    m_concentrationDecayRate(p->concentrationDecayRate)
{



}

void FlipSmokeSolver::applyBodyForces()
{
    float alpha = 0.01f;
    float beta = 1.f/m_ambientTemperature;
    for (int i = 0; i < m_sizeI + 1; i++)
    {
        for (int j = 0; j < m_sizeJ + 1; j++)
        {
            if(m_fluidVelocityGrid.velocityGridU().inBounds(i,j))
            {
                float tempAt = m_temperature.interpolateAt(static_cast<float>(i)-0.5f,j);
                float tempDiff = tempAt - m_ambientTemperature;
                float concentration = m_smokeConcentration.interpolateAt(static_cast<float>(i)-0.5f,j);
                float accelerationU = (alpha * concentration - beta * tempDiff) *
                                        m_globalAcceleration.x() * m_stepDt;
                m_fluidVelocityGrid.u(i,j) += accelerationU;
            }
            if(m_fluidVelocityGrid.velocityGridV().inBounds(i,j))
            {
                float tempAt = m_temperature.interpolateAt(i,static_cast<float>(j)-0.5f);
                float tempDiff = tempAt - m_ambientTemperature;
                float concentration = m_smokeConcentration.interpolateAt(i,static_cast<float>(j)-0.5f);
                float accelerationV = (alpha * concentration - beta * tempDiff) *
                                        m_globalAcceleration.y() * m_stepDt;
                m_fluidVelocityGrid.v(i,j) += accelerationV;
            }
        }
    }
}

void FlipSmokeSolver::centeredParamsToGrid()
{
    Grid2d<float> centeredWeights(m_sizeI,m_sizeJ,1e-10f);

    m_viscosityGrid.fill(0.f);
    m_temperature.fill(0.f);
    m_smokeConcentration.fill(0.f);
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
                        m_viscosityGrid.at(iIdx,jIdx) += weightCentered * p.viscosity;
                        m_temperature.at(iIdx,jIdx) += weightCentered * p.temperature;
                        m_smokeConcentration.at(iIdx,jIdx) += weightCentered * p.smokeConcentrartion;
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
                    m_viscosityGrid.at(i,j) /= centeredWeights.at(i,j);
                    m_temperature.at(i,j) /= centeredWeights.at(i,j);
                    m_temperature.at(i,j) = m_ambientTemperature +
                                    (m_temperature.at(i,j) - m_ambientTemperature) *
                                    std::exp(-m_temperatureDecayRate * m_stepDt);
                    m_smokeConcentration.at(i,j) /= centeredWeights.at(i,j);
                    m_smokeConcentration.at(i,j) *= std::exp(-m_concentrationDecayRate
                                                                   * m_stepDt);
                }
                else
                {
                    m_temperature.at(i,j) = m_ambientTemperature;
                }
            }
        }
    }
}

void FlipSmokeSolver::calcPressureRhs(std::vector<double> &rhs)
{
    double scale = 1/m_dx;

    for (int i = 0; i < m_sizeI; i++)
    {
        for (int j = 0; j < m_sizeJ; j++)
        {
            if (!m_materialGrid.isSolid(i,j))
            {
                rhs[linearIndex(i,j)] = -scale * divergenceAt(i,j);

                if(m_materialGrid.isSolid(i-1,j))
                {
                    rhs[linearIndex(i,j)] -= scale * static_cast<double>(m_fluidVelocityGrid.u(i,j) - 0);
                }
                if(m_materialGrid.isSolid(i+1,j))
                {
                    rhs[linearIndex(i,j)] += scale * static_cast<double>(m_fluidVelocityGrid.u(i+1,j) - 0);
                }

                if(m_materialGrid.isSolid(i,j-1))
                {
                    rhs[linearIndex(i,j)] -= scale * static_cast<double>(m_fluidVelocityGrid.v(i,j) - 0);
                }
                if(m_materialGrid.isSolid(i,j+1))
                {
                    rhs[linearIndex(i,j)] += scale * static_cast<double>(m_fluidVelocityGrid.v(i,j+1) - 0);
                }
            }
        }
    }
}

void FlipSmokeSolver::particleUpdate()
{
    FlipSolver::particleUpdate();
    for(int i = m_markerParticles.size() - 1; i >= 0; i--)
    {
        MarkerParticle &p = m_markerParticles[i];
        p.temperature = m_temperature.interpolateAt(p.position);
//        p.temperature = m_ambientTemperature +
//                (p.temperature - m_ambientTemperature) *
//                std::exp(-m_temperatureDecayRate * m_stepDt);
        p.smokeConcentrartion = m_smokeConcentration.interpolateAt(p.position);
//        p.smokeConcentrartion *= std::exp(-m_concentrationDecayRate * m_stepDt);
        p.testValue = p.smokeConcentrartion;
    }
}

void FlipSmokeSolver::afterTransfer()
{
    FlipSolver::afterTransfer();
    m_divergenceControl.fill(0.f);
    for (int i = 0; i < m_sizeI; i++)
    {
        for (int j = 0; j < m_sizeJ; j++)
        {
            if(m_materialGrid.isSource(i,j))
            {
                int emitterId = m_emitterId.at(i,j);
                m_smokeConcentration.at(i,j) = m_sources[emitterId].concentrartion();
                m_temperature.at(i,j) = m_sources[emitterId].temperature();
            }
        }
    }
}

void FlipSmokeSolver::reseedParticles()
{
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
                    addMarkerParticle(MarkerParticle{pos,velocity,viscosity,temp,conc});
                }
            }
        }
    }
}

void FlipSmokeSolver::applyPressuresToVelocityField(std::vector<double> &pressures)
{
    double scale = m_stepDt / (m_fluidDensity * m_dx);

    for (int i = m_sizeI - 1; i >= 0; i--)
    {
        for (int j = m_sizeJ - 1; j >= 0; j--)
        {
            int fluidIndex = linearIndex(i,j);
            int fluidIndexIM1 = linearIndex(i-1,j);
            int fluidIndexJM1 = linearIndex(i,j-1);
            double pCurrent = fluidIndex == -1 ? 0.0 : pressures[fluidIndex];
            //U part
            if(!m_materialGrid.isSolid(i-1,j) || !m_materialGrid.isSolid(i,j))
            {
                if(m_materialGrid.isSolid(i-1,j) || m_materialGrid.isSolid(i,j))
                {
                    m_fluidVelocityGrid.setU(i,j,0);//Solids are stationary
                }
                else
                {
                    double pIm1 = fluidIndexIM1 == -1 ? 0.0 : pressures[fluidIndexIM1];
                    m_fluidVelocityGrid.u(i,j) -= scale * (pCurrent - pIm1);
                }
            }
            else
            {
                m_fluidVelocityGrid.setUValidity(i,j,false);
            }

            //V part
            if(!m_materialGrid.isSolid(i,j-1) || !m_materialGrid.isSolid(i,j))
            {
                if(m_materialGrid.isSolid(i,j-1) || m_materialGrid.isSolid(i,j))
                {
                    m_fluidVelocityGrid.setV(i,j,0);//Solids are stationary
                }
                else
                {
                    double pJm1 = fluidIndexJM1 == -1 ? 0.0 : pressures[fluidIndexJM1];
                    m_fluidVelocityGrid.v(i,j) -= scale * (pCurrent - pJm1);
                }
            }
            else
            {
                m_fluidVelocityGrid.setVValidity(i,j,false);
            }
        }
    }

    if(anyNanInf(m_fluidVelocityGrid.velocityGridU().data()))
    {
        std::cout << "NaN or inf in U vector!\n" << std::flush;
    }

    if(anyNanInf(m_fluidVelocityGrid.velocityGridV().data()))
    {
        std::cout << "NaN or inf in V vector!\n" << std::flush;
    }
}

DynamicUpperTriangularSparseMatrix FlipSmokeSolver::getPressureProjectionMatrix()
{
    DynamicUpperTriangularSparseMatrix output = DynamicUpperTriangularSparseMatrix(cellCount(),7);

    double scale = m_stepDt / (m_fluidDensity * m_dx * m_dx);

    LinearIndexable2d& indexer = *dynamic_cast<LinearIndexable2d*>(this);

    for(int i = 0; i <  m_sizeI; i++)
    {
        for(int j = 0; j <  m_sizeJ; j++)
        {
            if(!m_materialGrid.isSolid(i,j))
            {
                //X Neighbors
                if(m_materialGrid.isFluid(i-1,j))
                {
                    output.addToAdiag(i,j,scale,  indexer);
                }else if(m_materialGrid.isEmpty(i-1,j))
                {
                    output.addToAdiag(i,j,scale,  indexer);
                }

                if(m_materialGrid.isFluid(i+1,j))
                {
                    output.addToAdiag(i,j,scale,  indexer);
                    output.setAx(i,j,-scale,  indexer);
                } else if(m_materialGrid.isEmpty(i+1,j))
                {
                    output.addToAdiag(i,j,scale,  indexer);
                }

                //Y Neighbors
                if(m_materialGrid.isFluid(i,j-1))
                {
                    output.addToAdiag(i,j,scale,  indexer);
                }else if(m_materialGrid.isEmpty(i,j-1))
                {
                    output.addToAdiag(i,j,scale,  indexer);
                }

                if(m_materialGrid.isFluid(i,j+1))
                {
                    output.addToAdiag(i,j,scale,  indexer);
                    output.setAy(i,j,-scale,  indexer);
                } else if(m_materialGrid.isEmpty(i,j+1))
                {
                    output.addToAdiag(i,j,scale,  indexer);
                }
            }
        }
    }

    return output;
}

const Grid2d<float> FlipSmokeSolver::smokeConcentration() const
{
    return m_smokeConcentration;
}

const Grid2d<float> FlipSmokeSolver::temperature() const
{
    return m_temperature;
}
