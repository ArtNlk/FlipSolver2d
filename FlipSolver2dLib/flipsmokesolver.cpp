#include "flipsmokesolver.h"

#include "mathfuncs.h"
#include "simsettings.h"

FlipSmokeSolver::FlipSmokeSolver(int extrapRadius, bool vonNeumannNeighbors):
    FlipSolver(extrapRadius, vonNeumannNeighbors)
{
}

void FlipSmokeSolver::step()
{
    updateSdf();
    Grid2d<int> particleCounts(m_grid.sizeI(), m_grid.sizeJ());
    m_grid.updateLinearFluidViscosityMapping();
    countParticles();
    reseedParticles();
    updateMaterials();
    particleToGrid();
    extrapolateVelocityField(m_grid.fluidVelocityGridU(),m_grid.knownFluidFlagsGridU(),10);
    extrapolateVelocityField(m_grid.fluidVelocityGridV(),m_grid.knownFluidFlagsGridV(),10);
    m_grid.savedFluidVelocityGrid().velocityGridU() = m_grid.fluidVelocityGridU();
    m_grid.savedFluidVelocityGrid().velocityGridV() = m_grid.fluidVelocityGridV();
    applyBodyForces();
    project();
    updateVelocityFromSolids();
//    applyViscosity();
//    project();
    extrapolateVelocityField(m_grid.fluidVelocityGridU(),m_grid.knownFluidFlagsGridU(),10);
    extrapolateVelocityField(m_grid.fluidVelocityGridV(),m_grid.knownFluidFlagsGridV(),10);
    particleUpdate();
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
            if(m_grid.fluidVelocityGridU().inBounds(i,j))
            {
                float tempAt = m_grid.temperatureAt(i,static_cast<float>(j) + 0.5f);
                float tempDiff = tempAt - SimSettings::ambientTemp();
                float concentration = m_grid.smokeConcentrationAt(i,static_cast<float>(j) + 0.5f);
                float accelerationU = (alpha * concentration - beta * tempDiff) *
                                        SimSettings::globalAcceleration().x() * SimSettings::stepDt();
                m_grid.fluidVelocityGridU().at(i,j) += accelerationU;
            }
            if(m_grid.fluidVelocityGridV().inBounds(i,j))
            {
                float tempAt = m_grid.temperatureAt(static_cast<float>(i) + 0.5f,j);
                float tempDiff = tempAt - SimSettings::ambientTemp();
                float concentration = m_grid.smokeConcentrationAt(static_cast<float>(i) + 0.5f,j);
                float accelerationV = (alpha * concentration - beta * tempDiff) *
                                        SimSettings::globalAcceleration().y() * SimSettings::stepDt();
                m_grid.fluidVelocityGridV().at(i,j) += accelerationV;
            }
        }
    }
}

void FlipSmokeSolver::particleToGrid()
{
    ASSERT(m_grid.fluidVelocityGridU().sizeI() == m_grid.fluidVelocityGridU().sizeI() && m_grid.fluidVelocityGridU().sizeJ() == m_grid.fluidVelocityGridU().sizeJ());
    ASSERT(m_grid.fluidVelocityGridV().sizeI() == m_grid.fluidVelocityGridV().sizeI() && m_grid.fluidVelocityGridV().sizeJ() == m_grid.fluidVelocityGridV().sizeJ());

    Grid2d<float> uWeights(m_grid.fluidVelocityGridU().sizeI(),m_grid.fluidVelocityGridU().sizeJ(),1e-10f);
    Grid2d<float> vWeights(m_grid.fluidVelocityGridV().sizeI(),m_grid.fluidVelocityGridV().sizeJ(),1e-10f);
    Grid2d<float> centeredWeights(m_grid.sizeI(),m_grid.sizeJ(),1e-10f);

    m_grid.fluidVelocityGridU().fill(0.f);
    m_grid.fluidVelocityGridV().fill(0.f);
    m_grid.viscosityGrid().fill(0.f);
    m_grid.temperatureGrid().fill(0.f);
    m_grid.smokeConcentrationGrid().fill(0.f);
    m_grid.knownFluidFlagsGridU().fill(false);
    m_grid.knownFluidFlagsGridV().fill(false);
    m_grid.knownFlagsCenteredParams().fill(false);

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
                float weightCentered = simmath::quadraticBSpline(p.position.x() - (iIdx),
                                                     p.position.y() - (jIdx));
                if(uWeights.inBounds(iIdx,jIdx))
                {
                    if(std::abs(weightU) > 1e-9f && std::abs(weightV) > 1e-9f)
                    {
                        uWeights.at(iIdx,jIdx) += weightU;
                        m_grid.fluidVelocityGridU().at(iIdx,jIdx) += weightU * (p.velocity.x() * SimSettings::dx());
                        m_grid.knownFluidFlagsGridU().at(iIdx,jIdx) = true;
                    }
                }

                if(vWeights.inBounds(iIdx,jIdx))
                {
                    if(std::abs(weightU) > 1e-9f && std::abs(weightV) > 1e-9f)
                    {
                        vWeights.at(iIdx,jIdx) += weightV;
                        m_grid.fluidVelocityGridV().at(iIdx,jIdx) += weightV * (p.velocity.y() * SimSettings::dx());
                        m_grid.knownFluidFlagsGridV().at(iIdx,jIdx) = true;
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
            if(m_grid.fluidVelocityGridU().inBounds(i,j))
            {
                m_grid.fluidVelocityGridU().at(i,j) /= uWeights.at(i,j);
            }
            if(m_grid.fluidVelocityGridV().inBounds(i,j))
            {
                m_grid.fluidVelocityGridV().at(i,j) /= vWeights.at(i,j);
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
                rhs[m_grid.linearIndex(i,j)] = -scale * static_cast<double>(m_grid.getFluidU(i+1,j) - m_grid.getFluidU(i,j)
                                                              +m_grid.getFluidV(i,j+1) - m_grid.getFluidV(i,j))
                                                                + m_grid.divergenceControl(i,j);

                if(m_grid.isSolid(i-1,j))
                {
                    rhs[m_grid.linearIndex(i,j)] -= scale * static_cast<double>(m_grid.getFluidU(i,j) - 0);
                }
                if(m_grid.isSolid(i+1,j))
                {
                    rhs[m_grid.linearIndex(i,j)] += scale * static_cast<double>(m_grid.getFluidU(i+1,j) - 0);
                }

                if(m_grid.isSolid(i,j-1))
                {
                    rhs[m_grid.linearIndex(i,j)] -= scale * static_cast<double>(m_grid.getFluidV(i,j) - 0);
                }
                if(m_grid.isSolid(i,j+1))
                {
                    rhs[m_grid.linearIndex(i,j)] += scale * static_cast<double>(m_grid.getFluidV(i,j+1) - 0);
                }
            }
        }
    }
}

void FlipSmokeSolver::applyPressuresToVelocityField(std::vector<double> &pressures)
{
    double scale = SimSettings::stepDt() / (SimSettings::fluidDensity() * SimSettings::dx());

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
                    m_grid.setFluidU(i,j,0);//Solids are stationary
                }
                else
                {
                    m_grid.fluidVelocityGridU().at(i,j) -= scale * (fluidIndex == -1 ? 0.0 : pressures[fluidIndex] - (fluidIndexIM1 == -1 ? 0.0 : pressures[fluidIndexIM1]));
                }
                m_grid.knownFluidFlagsGridU().at(i,j) = true;
            }
            else
            {
                m_grid.knownFluidFlagsGridU().at(i,j) = false;
            }

            //V part
            if(!m_grid.isSolid(i,j-1) || !m_grid.isSolid(i,j))
            {
                if(m_grid.isSolid(i,j-1) || m_grid.isSolid(i,j))
                {
                    m_grid.fluidVelocityGridV().at(i,j) = 0;//Solids are stationary
                }
                else
                {
                    m_grid.fluidVelocityGridV().at(i,j) -= scale * (fluidIndex == -1 ? 0.0 : pressures[fluidIndex] - (fluidIndexJM1 == -1 ? 0.0 : pressures[fluidIndexJM1]));
                }
                m_grid.knownFluidFlagsGridV().at(i,j) = true;
            }
            else
            {
                m_grid.knownFluidFlagsGridV().at(i,j) = false;
            }
        }
    }

    if(anyNanInf(m_grid.fluidVelocityGridU().data()))
    {
        std::cout << "NaN or inf in U vector!\n" << std::flush;
    }

    if(anyNanInf(m_grid.fluidVelocityGridV().data()))
    {
        std::cout << "NaN or inf in V vector!\n" << std::flush;
    }
}

DynamicUpperTriangularSparseMatrix FlipSmokeSolver::getPressureProjectionMatrix()
{
    DynamicUpperTriangularSparseMatrix output = DynamicUpperTriangularSparseMatrix(m_grid.cellCount(),7);

    double scale = SimSettings::stepDt() / (SimSettings::fluidDensity() * SimSettings::dx() * SimSettings::dx());

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
