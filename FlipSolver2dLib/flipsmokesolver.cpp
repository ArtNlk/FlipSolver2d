#include "flipsmokesolver.h"

#include "flipsolver2d.h"
#include "mathfuncs.h"
#include "simsettings.h"

FlipSmokeSolver::FlipSmokeSolver(int sizeI, int sizeJ, int extrapNeighborRadius, bool vonNeumannNeighbors):
    FlipSolver(sizeI, sizeJ, extrapNeighborRadius, vonNeumannNeighbors),
    m_temperature(sizeI, sizeJ, SimSettings::ambientTemp(), OOBStrategy::OOB_CONST, SimSettings::ambientTemp()),
    m_smokeConcentration(sizeI, sizeJ, 0.f, OOBStrategy::OOB_CONST, 0.f)
{
}

void FlipSmokeSolver::applyBodyForces()
{
    float alpha = 0.01f;
    float beta = 1.f/SimSettings::ambientTemp();
    for (int i = 0; i < m_sizeI + 1; i++)
    {
        for (int j = 0; j < m_sizeJ + 1; j++)
        {
            if(m_fluidVelocityGrid.velocityGridU().inBounds(i,j))
            {
                float tempAt = m_temperature.interpolateAt(i,static_cast<float>(j) + 0.5f);
                float tempDiff = tempAt - SimSettings::ambientTemp();
                float concentration = m_smokeConcentration.interpolateAt(i,static_cast<float>(j) + 0.5f);
                float accelerationU = (alpha * concentration - beta * tempDiff) *
                                        SimSettings::globalAcceleration().x() * SimSettings::stepDt();
                m_fluidVelocityGrid.u(i,j) += accelerationU;
            }
            if(m_fluidVelocityGrid.velocityGridV().inBounds(i,j))
            {
                float tempAt = m_temperature.interpolateAt(static_cast<float>(i) + 0.5f,j);
                float tempDiff = tempAt - SimSettings::ambientTemp();
                float concentration = m_smokeConcentration.interpolateAt(static_cast<float>(i) + 0.5f,j);
                float accelerationV = (alpha * concentration - beta * tempDiff) *
                                        SimSettings::globalAcceleration().y() * SimSettings::stepDt();
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
                if(inBounds(iIdx,jIdx))
                {
                    float weightCentered = simmath::quadraticBSpline(p.position.x() - (iIdx),
                                                         p.position.y() - (jIdx));
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
            if(centeredWeights.inBounds(i,j))
            {
                if(m_knownCenteredParams.at(i,j))
                {
                    m_viscosityGrid.at(i,j) /= centeredWeights.at(i,j);
                    m_temperature.at(i,j) /= centeredWeights.at(i,j);
                    m_smokeConcentration.at(i,j) /= centeredWeights.at(i,j);
                }
                else
                {
                    m_temperature.at(i,j) = SimSettings::ambientTemp();
                }
            }
        }
    }
}

void FlipSmokeSolver::calcPressureRhs(std::vector<double> &rhs)
{
    double scale = 1/SimSettings::dx();

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

void FlipSmokeSolver::applyPressuresToVelocityField(std::vector<double> &pressures)
{
    double scale = SimSettings::stepDt() / (SimSettings::fluidDensity() * SimSettings::dx());

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

    double scale = SimSettings::stepDt() / (SimSettings::fluidDensity() * SimSettings::dx() * SimSettings::dx());

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
