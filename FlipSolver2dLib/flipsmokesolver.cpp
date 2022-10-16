#include "flipsmokesolver.h"

#include "mathfuncs.h"
#include "simsettings.h"

FlipSmokeSolver::FlipSmokeSolver(int extrapRadius, bool vonNeumannNeighbors):
    FlipSolver(extrapRadius, vonNeumannNeighbors)
{
    //m_grid.temperatureGrid().fill(SimSettings::ambientTemp());
}

void FlipSmokeSolver::step()
{
    Grid2d<int> particleCounts(m_grid.sizeI(), m_grid.sizeJ());
    m_grid.updateLinearFluidViscosityMapping();
    countParticles(particleCounts);
    reseedParticles(particleCounts);
    updateMaterialsFromParticles(particleCounts);
    particleToGrid();
    extrapolateVelocityField(1);
    Grid2d<float> prevU = m_grid.velocityGridU();
    Grid2d<float> prevV = m_grid.velocityGridV();
    applyBodyForces();
    project();
    updateVelocityFromSolids();
//    applyViscosity();
//    project();
    extrapolateVelocityField(1);
    particleUpdate(prevU, prevV);
    advect();
}

void FlipSmokeSolver::applyBodyForces()
{
    float alpha = 0.01f;
    float beta = 1.f/SimSettings::ambientTemp();
    for (int i = 0; i < m_grid.sizeI() + 1; i++)
    {
        for (int j = 0; j < m_grid.sizeJ() + 1; j++)
        {
            if(m_grid.velocityGridU().inBounds(i,j))
            {
                float tempAt = m_grid.temperatureAt(i,static_cast<float>(j) + 0.5f);
                float tempDiff = tempAt - SimSettings::ambientTemp();
                float concentration = m_grid.smokeConcentrationAt(i,static_cast<float>(j) + 0.5f);
                float accelerationU = (alpha * concentration - beta * tempDiff) *
                                        SimSettings::globalAcceleration().x() * SimSettings::stepDt();
                m_grid.velocityGridU().at(i,j) += accelerationU;
            }
            if(m_grid.velocityGridV().inBounds(i,j))
            {
                float tempAt = m_grid.temperatureAt(static_cast<float>(i) + 0.5f,j);
                float tempDiff = tempAt - SimSettings::ambientTemp();
                float concentration = m_grid.smokeConcentrationAt(static_cast<float>(i) + 0.5f,j);
                float accelerationV = (alpha * concentration - beta * tempDiff) *
                                        SimSettings::globalAcceleration().y() * SimSettings::stepDt();
                m_grid.velocityGridV().at(i,j) += accelerationV;
            }
        }
    }
}

void FlipSmokeSolver::particleToGrid()
{
    ASSERT(m_grid.velocityGridU().sizeI() == m_grid.velocityGridU().sizeI() && m_grid.velocityGridU().sizeJ() == m_grid.velocityGridU().sizeJ());
    ASSERT(m_grid.velocityGridV().sizeI() == m_grid.velocityGridV().sizeI() && m_grid.velocityGridV().sizeJ() == m_grid.velocityGridV().sizeJ());

    Grid2d<float> uWeights(m_grid.velocityGridU().sizeI(),m_grid.velocityGridU().sizeJ(),1e-10f);
    Grid2d<float> vWeights(m_grid.velocityGridV().sizeI(),m_grid.velocityGridV().sizeJ(),1e-10f);
    Grid2d<float> centeredWeights(m_grid.sizeI(),m_grid.sizeJ(),1e-10f);

    m_grid.velocityGridU().fill(0.f);
    m_grid.velocityGridV().fill(0.f);
    m_grid.viscosityGrid().fill(0.f);
    m_grid.temperatureGrid().fill(0.f);
    m_grid.smokeConcentrationGrid().fill(0.f);
    m_grid.knownFlagsGridU().fill(false);
    m_grid.knownFlagsGridV().fill(false);
    m_grid.knownFlagsCenteredParams().fill(false);

    for(MarkerParticle &p : m_markerParticles)
    {
        int i = math::integr(p.position.x());
        int j = math::integr(p.position.y());
        //Run over all cells that this particle might affect
        for (int iOffset = -3; iOffset < 3; iOffset++)
        {
            for (int jOffset = -3; jOffset < 3; jOffset++)
            {
                int iIdx = i+iOffset;
                int jIdx = j+jOffset;
                float weightU = math::quadraticBSpline(p.position.x() - (iIdx),
                                                     p.position.y() - (static_cast<float>(jIdx) + 0.5f));

                float weightV = math::quadraticBSpline(p.position.x() - (static_cast<float>(iIdx) + 0.5f),
                                                     p.position.y() - (jIdx));
                float weightCentered = math::quadraticBSpline(p.position.x() - (iIdx),
                                                     p.position.y() - (jIdx));
                if(uWeights.inBounds(iIdx,jIdx))
                {
                    if(std::abs(weightU) > 1e-9f && std::abs(weightV) > 1e-9f)
                    {
                        uWeights.at(iIdx,jIdx) += weightU;
                        m_grid.velocityGridU().at(iIdx,jIdx) += weightU * (p.velocity.x() * SimSettings::dx());
                        m_grid.knownFlagsGridU().at(iIdx,jIdx) = true;
                    }
                }

                if(vWeights.inBounds(iIdx,jIdx))
                {
                    if(std::abs(weightU) > 1e-9f && std::abs(weightV) > 1e-9f)
                    {
                        vWeights.at(iIdx,jIdx) += weightV;
                        m_grid.velocityGridV().at(iIdx,jIdx) += weightV * (p.velocity.y() * SimSettings::dx());
                        m_grid.knownFlagsGridV().at(iIdx,jIdx) = true;
                    }
                }

                if(centeredWeights.inBounds(iIdx,jIdx))
                {
                    if(std::abs(weightCentered) > 1e-6f)
                    {
                        centeredWeights.at(iIdx,jIdx) += weightCentered;
                        m_grid.viscosityGrid().at(iIdx,jIdx) += weightCentered * p.viscosity;
                        m_grid.temperatureGrid().at(iIdx,jIdx) += weightCentered * p.temperature;
                        m_grid.smokeConcentrationGrid().at(iIdx,jIdx) += weightCentered * p.smokeConcentrartion;
                        m_grid.knownFlagsCenteredParams().at(iIdx,jIdx) = true;
                    }
                }
            }
        }
    }

    for (int i = 0; i < m_grid.sizeI() + 1; i++)
    {
        for (int j = 0; j < m_grid.sizeJ() + 1; j++)
        {
            if(m_grid.velocityGridU().inBounds(i,j))
            {
                m_grid.velocityGridU().at(i,j) /= uWeights.at(i,j);
            }
            if(m_grid.velocityGridV().inBounds(i,j))
            {
                m_grid.velocityGridV().at(i,j) /= vWeights.at(i,j);
            }
            if(centeredWeights.inBounds(i,j))
            {
                if(m_grid.knownFlagsCenteredParams().at(i,j))
                {
                    m_grid.viscosityGrid().at(i,j) /= centeredWeights.at(i,j);
                    m_grid.temperatureGrid().at(i,j) /= centeredWeights.at(i,j);
                    m_grid.smokeConcentrationGrid().at(i,j) /= centeredWeights.at(i,j);
                }
                else
                {
                    m_grid.temperatureGrid().at(i,j) = SimSettings::ambientTemp();
                }
            }
        }
    }
}

void FlipSmokeSolver::project()
{
    m_grid.updateLinearFluidViscosityMapping();
    std::vector<double> rhs(m_grid.cellCount(),0.0);
    calcPressureRhs(rhs);
    //debug() << "Calculated rhs: " << rhs;
    DynamicUpperTriangularSparseMatrix mat = getPressureProjectionMatrix();
    std::vector<double> pressures(m_grid.cellCount(),0.0);
    if(!m_pcgSolver.solve(mat,pressures,rhs,200))
    {
        std::cout << "PCG Solver pressure: Iteration limit exhaustion!\n";
    }

    if(anyNanInf(pressures))
    {
        std::cout << "NaN or inf in pressures vector!\n" << std::flush;
    }

    //debug() << "pressures = " << pressures;

    applyPressuresToVelocityField(pressures);
}

void FlipSmokeSolver::calcPressureRhs(std::vector<double> &rhs)
{
    double scale = 1/SimSettings::dx();

    for (int i = 0; i < m_grid.sizeI(); i++)
    {
        for (int j = 0; j < m_grid.sizeJ(); j++)
        {
            if (!m_grid.isSolid(i,j))
            {
                rhs[m_grid.linearIndex(i,j)] = -scale * static_cast<double>(m_grid.getU(i+1,j) - m_grid.getU(i,j)
                                                              +m_grid.getV(i,j+1) - m_grid.getV(i,j))
                                                                + m_grid.divergenceControl(i,j);

                if(m_grid.isSolid(i-1,j))
                {
                    rhs[m_grid.linearIndex(i,j)] -= scale * static_cast<double>(m_grid.getU(i,j) - 0);
                }
                if(m_grid.isSolid(i+1,j))
                {
                    rhs[m_grid.linearIndex(i,j)] += scale * static_cast<double>(m_grid.getU(i+1,j) - 0);
                }

                if(m_grid.isSolid(i,j-1))
                {
                    rhs[m_grid.linearIndex(i,j)] -= scale * static_cast<double>(m_grid.getV(i,j) - 0);
                }
                if(m_grid.isSolid(i,j+1))
                {
                    rhs[m_grid.linearIndex(i,j)] += scale * static_cast<double>(m_grid.getV(i,j+1) - 0);
                }
            }
        }
    }
}

void FlipSmokeSolver::applyPressuresToVelocityField(std::vector<double> &pressures)
{
    double scale = SimSettings::stepDt() / (SimSettings::density() * SimSettings::dx());

    for (int i = m_grid.sizeI() - 1; i >= 0; i--)
    {
        for (int j = m_grid.sizeJ() - 1; j >= 0; j--)
        {
            int fluidIndex = m_grid.linearIndex(i,j);
            int fluidIndexIM1 = m_grid.linearIndex(i-1,j);
            int fluidIndexJM1 = m_grid.linearIndex(i,j-1);
            //U part
            if(!m_grid.isSolid(i-1,j) || !m_grid.isSolid(i,j))
            {
                if(m_grid.isSolid(i-1,j) || m_grid.isSolid(i,j))
                {
                    m_grid.setU(i,j,0);//Solids are stationary
                }
                else
                {
                    m_grid.velocityGridU().at(i,j) -= scale * (fluidIndex == -1 ? 0.0 : pressures[fluidIndex] - (fluidIndexIM1 == -1 ? 0.0 : pressures[fluidIndexIM1]));
                }
                m_grid.knownFlagsGridU().at(i,j) = true;
            }
            else
            {
                m_grid.knownFlagsGridU().at(i,j) = false;
            }

            //V part
            if(!m_grid.isSolid(i,j-1) || !m_grid.isSolid(i,j))
            {
                if(m_grid.isSolid(i,j-1) || m_grid.isSolid(i,j))
                {
                    m_grid.velocityGridV().at(i,j) = 0;//Solids are stationary
                }
                else
                {
                    m_grid.velocityGridV().at(i,j) -= scale * (fluidIndex == -1 ? 0.0 : pressures[fluidIndex] - (fluidIndexJM1 == -1 ? 0.0 : pressures[fluidIndexJM1]));
                }
                m_grid.knownFlagsGridV().at(i,j) = true;
            }
            else
            {
                m_grid.knownFlagsGridV().at(i,j) = false;
            }
        }
    }

    if(anyNanInf(m_grid.velocityGridU().data()))
    {
        std::cout << "NaN or inf in U vector!\n" << std::flush;
    }

    if(anyNanInf(m_grid.velocityGridV().data()))
    {
        std::cout << "NaN or inf in V vector!\n" << std::flush;
    }
}

DynamicUpperTriangularSparseMatrix FlipSmokeSolver::getPressureProjectionMatrix()
{
    DynamicUpperTriangularSparseMatrix output = DynamicUpperTriangularSparseMatrix(m_grid.cellCount(),7);

    double scale = SimSettings::stepDt() / (SimSettings::density() * SimSettings::dx() * SimSettings::dx());

    for(int i = 0; i <  m_grid.sizeI(); i++)
    {
        for(int j = 0; j <  m_grid.sizeJ(); j++)
        {
            if(!m_grid.isSolid(i,j))
            {
                //X Neighbors
                if( m_grid.isFluid(i-1,j))
                {
                    output.addToAdiag(i,j,scale,  m_grid);
                }else if( m_grid.isEmpty(i-1,j))
                {
                    output.addToAdiag(i,j,scale,  m_grid);
                }

                if( m_grid.isFluid(i+1,j))
                {
                    output.addToAdiag(i,j,scale,  m_grid);
                    output.setAx(i,j,-scale,  m_grid);
                } else if( m_grid.isEmpty(i+1,j))
                {
                    output.addToAdiag(i,j,scale,  m_grid);
                }

                //Y Neighbors
                if( m_grid.isFluid(i,j-1))
                {
                    output.addToAdiag(i,j,scale,  m_grid);
                }else if( m_grid.isEmpty(i,j-1))
                {
                    output.addToAdiag(i,j,scale,  m_grid);
                }

                if( m_grid.isFluid(i,j+1))
                {
                    output.addToAdiag(i,j,scale,  m_grid);
                    output.setAy(i,j,-scale,  m_grid);
                } else if( m_grid.isEmpty(i,j+1))
                {
                    output.addToAdiag(i,j,scale,  m_grid);
                }
            }
        }
    }

    return output;
}
