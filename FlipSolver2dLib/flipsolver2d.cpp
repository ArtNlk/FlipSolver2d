#include "flipsolver2d.h"

#include "simsettings.h"

FlipSolver::FlipSolver(int sizeX, int sizeY, double fluidDensity, double timestepSize,  double sideLength) :
    m_grid(sizeX, sizeY)
{
    SimSettings::density() = fluidDensity;
    SimSettings::dt() = timestepSize;
    SimSettings::dx() = sideLength;
}

void FlipSolver::project()
{
    std::vector<double> rhs(0,m_grid.sizeI() * m_grid.sizeJ());
    calcRhs(rhs);
    SparseMatrix mat = DynamicSparseMatrix(DynamicSparseMatrix(m_grid));
    std::vector<double> pressures(0,m_grid.sizeI() * m_grid.sizeJ());
    if(!m_pcgSolver.solve(mat,m_grid,pressures,rhs))
    {
        //TODO report iteration limit exhaustion
    }

    double scale = SimSettings::dt() / (SimSettings::density() * SimSettings::dx());

    for (int i = m_grid.sizeI(); i >= 1; i--)
    {
        for (int j = m_grid.sizeJ(); j >= 1; j--)
        {
            //U part
            if(m_grid.getMaterial(i-1,j) == FluidCellMaterial::FLUID || m_grid.getMaterial(i,j) == FluidCellMaterial::FLUID)
            {
                if(m_grid.getMaterial(i-1,j) == FluidCellMaterial::SOLID || m_grid.getMaterial(i,j) == FluidCellMaterial::SOLID)
                {
                    m_grid.at(i,j).setU(0);//Solids are stationary
                }
                else
                {
                    m_grid.at(i,j).U() -= scale * (pressures[m_grid.linearIndex(i,j)] - pressures[m_grid.linearIndex(i-1,j)]);
                }
            }
            else
            {
                m_grid.at(i,j).setUnknownU();
            }

            //V part
            if(m_grid.getMaterial(i,j-1) == FluidCellMaterial::FLUID || m_grid.getMaterial(i,j) == FluidCellMaterial::FLUID)
            {
                if(m_grid.getMaterial(i,j-1) == FluidCellMaterial::SOLID || m_grid.getMaterial(i,j) == FluidCellMaterial::SOLID)
                {
                    m_grid.at(i,j).setV(0);//Solids are stationary
                }
                else
                {
                    m_grid.at(i,j).V() -= scale * (pressures[m_grid.linearIndex(i,j)] - pressures[m_grid.linearIndex(i,j-1)]);
                }
            }
            else
            {
                m_grid.at(i,j).setUnknownV();
            }
        }
    }
}

void FlipSolver::calcRhs(std::vector<double> &rhs)
{
    double scale = 1/SimSettings::dx();

    for (int i = 1; i < m_grid.sizeI(); i++)
    {
        for (int j = 1; j < m_grid.sizeJ(); j++)
        {
            if (m_grid.getMaterial(i,j) == FluidCellMaterial::FLUID)
            {
                rhs[m_grid.linearIndex(i,j)] = -scale * (m_grid.getU(i+1,j) - m_grid.getU(i,j)
                                                          +m_grid.getV(i,j+1) - m_grid.getV(i,j));
            }
        }
    }
}
