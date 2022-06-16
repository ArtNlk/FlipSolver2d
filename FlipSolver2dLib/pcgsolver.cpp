#include "pcgsolver.h"

#include <algorithm>

#include <sstream>

#include "vmath.h"
#include "cmath"
#include "simsettings.h"
#include "logger.h"

const double PCGSolver::m_tol = 1.0e-14;

PCGSolver::PCGSolver() :
    m_precondCache()
{

}

bool PCGSolver::solve(const UpperTriangularMatrix &matrix, MACFluidGrid &grid, std::vector<double> &result, const std::vector<double> &vec, int iterLimit)
{
    result.assign(result.size(),0);
    if (vmath::isZero(vec))
    {
        return false;
    }
    calcPrecond(matrix,grid);
    m_residual = vec;
    m_aux = vec;
    applyICPrecond(matrix,m_residual,m_aux,grid);
    m_search = m_aux;
    double sigma = vmath::dot(m_aux,m_residual);
    double err = 0.0;
    for (int i = 0; i < iterLimit; i++)
    {
        m_aux = matrix * m_search;
        double alpha = sigma/vmath::dot(m_aux,m_search);
        vmath::addMul(result,result,m_search,alpha);
        vmath::subMul(m_residual,m_residual,m_aux,alpha);
        err = vmath::maxAbs(m_residual);
//        if(i % 5 == 0)
//        {
//            std::cout << "Solver: " << i << " : " << err << "\n";
//            debug() << "Solver: " << i << " : " << err;
//        }
        if (err <= m_tol)
        {
            debug() << "Solver done, iter = " << i << " err = " << err;
            std::cout << "Solver done, iter = " << i << " err = " << err << '\n';
            return true;
        }
        m_aux = m_residual;
        applyICPrecond(matrix,m_residual,m_aux,grid);
        double newSigma = vmath::dot(m_aux,m_residual);
        double beta = newSigma/sigma;
        vmath::addMul(m_search,m_aux,m_search,beta);
        sigma = newSigma;
    }

    debug() << "Solver iter exhaustion, err = " << err;
    return false;
}

void PCGSolver::applyICPrecond(const UpperTriangularMatrix &matrix, std::vector<double> &in, std::vector<double> &out, MACFluidGrid &grid)
{
    double t = 0;
    std::vector<double> temp(in.size(),0);
    for (int i = 0; i < grid.sizeI(); i++)
    {
        for (int j = 0; j < grid.sizeJ(); j++)
        {
            if (grid.isFluid(i,j))
            {
                double v1 = grid.linearFluidIndex(i-1,j) == -1? 0.0 : matrix.Ax(i-1,j,grid) * precond(matrix,i-1,j,grid) * temp[grid.linearFluidIndex(i-1,j)];
                double v2 = grid.linearFluidIndex(i,j-1) == -1? 0.0 : matrix.Ay(i,j-1,grid) * precond(matrix,i,j-1,grid) * temp[grid.linearFluidIndex(i,j-1)];

                t = in[grid.linearFluidIndex(i,j)] - v1 - v2;
                temp[grid.linearFluidIndex(i,j)] = t*precond(matrix,i,j,grid);
            }
        }
    }

    for (int i = grid.sizeI() - 1; i >= 0; i--)
    {
        for (int j = grid.sizeJ() - 1; j >= 0; j--)
        {
            if (grid.isFluid(i,j))
            {
                double v1 = grid.linearFluidIndex(i+1,j) == -1? 0.0 : matrix.Ax(i,j,grid) * precond(matrix,i,j,grid) * out[grid.linearFluidIndex(i+1,j)];
                double v2 = grid.linearFluidIndex(i,j+1) == -1? 0.0 : matrix.Ay(i,j,grid) * precond(matrix,i,j,grid) * out[grid.linearFluidIndex(i,j+1)];

                t = temp[grid.linearFluidIndex(i,j)] - v1 - v2;
                out[grid.linearFluidIndex(i,j)] = t*precond(matrix,i,j,grid);
            }
        }
    }
}

void PCGSolver::calcPrecond(const UpperTriangularMatrix &matrix, MACFluidGrid &grid)
{
    m_precondCache.clear();
    for (int i = 0; i < grid.sizeI(); i++)
    {
        for (int j = 0; j < grid.sizeJ(); j++)
        {
            if (grid.isFluid(i,j))
            {
                precond(matrix,i,j,grid);
            }

        }
    }
}

double PCGSolver::precond(const UpperTriangularMatrix &m, int i, int j, MACFluidGrid& grid)
{
    if(i < 0 || j < 0)
    {
        return 0.0;
    }
    double temp;
    std::unordered_map<int,double>::iterator iter = m_precondCache.find(grid.linearFluidIndex(i,j));
    if(iter != m_precondCache.end())
    {
        return iter->second;
    }

    const double tuning = 0.97;
    const double safety = 0.25;

    temp = m.Adiag(i,j,grid)
            - (m.Ax(i-1,j,grid)*precond(m,i-1,j,grid))*(m.Ax(i-1,j,grid)*precond(m,i-1,j,grid))
            - (m.Ay(i,j-1,grid)*precond(m,i,j-1,grid))*(m.Ay(i,j-1,grid)*precond(m,i,j-1,grid))

            - tuning*(m.Ax(i-1,j,grid)*(m.Ay(i-1,j,grid))*precond(m,i-1,j,grid)*precond(m,i-1,j,grid)
                      + m.Ay(i,j-1,grid)*(m.Ax(i,j-1,grid))*precond(m,i,j-1,grid)*precond(m,i,j-1,grid));

    if(temp < safety*m.Adiag(i,j,grid))
    {
        temp = m.Adiag(i,j,grid);
    }

    if(std::abs(temp) > 10e-10)
    {
        temp = 1/std::sqrt(temp);
    }
    else
    {
        temp = 0.0;
    }
    m_precondCache[grid.linearFluidIndex(i,j)] = temp;

    return temp;
}
