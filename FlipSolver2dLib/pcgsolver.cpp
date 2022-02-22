#include "pcgsolver.h"

#include <algorithm>

#include "vmath.h"
#include "cmath"

const double PCGSolver::m_tol = 1.0e-6;

PCGSolver::PCGSolver() :
    m_precondCache(1,1)
{

}

bool PCGSolver::solve(const SparseMatrix &matrix, MACFluidGrid &grid, std::vector<double> &result, const std::vector<double> &vec, int iterLimit)
{
    result.assign(matrix.rowCount(),0);
    if (vmath::isZero(vec))
    {
        return false;
    }
    calcPrecond(matrix,grid);
    m_residual = vec;
    m_aux = m_residual;
    applyICPrecond(matrix,m_aux,grid);
    m_search = m_aux;
    double sigma = vmath::dot(m_aux,m_residual);

    for (int i = 0; i < iterLimit; i++)
    {
        m_aux = matrix * m_search;
        double alpha = sigma/vmath::dot(m_aux,m_search);
        vmath::addMul(result,result,m_search,alpha);
        vmath::subMul(m_residual,m_residual,m_aux,alpha);
        if (vmath::maxAbs(m_residual) <= m_tol)
        {
            return true;
        }
        m_aux = m_residual;
        applyICPrecond(matrix,m_aux,grid);
        double newSigma = vmath::dot(m_aux,m_residual);
        double beta = newSigma/sigma;
        vmath::addMul(m_search,m_aux,m_search,beta);
        sigma = newSigma;
    }

    return false;
}

void PCGSolver::applyICPrecond(const SparseMatrix &matrix, std::vector<double> &vector, MACFluidGrid &grid)
{
    double t = 0;
    std::vector<double> temp(vector.size(),0);
    for (int i = 1; i < grid.sizeI(); i++)
    {
        for (int j = 1; j < grid.sizeJ(); j++)
        {
            if (grid.getMaterial(i,j) == FluidCellMaterial::FLUID)
            {
                t = vector[matrix.linearIndex(i-1,j)] - (matrix.Ax(i-1,j) * precond(matrix,i-1,j) * temp[matrix.linearIndex(i-1,j)])
                                                        - (matrix.Ax(i,j-1) * precond(matrix,i,j-1) * temp[matrix.linearIndex(i,j-1)]);
                temp[matrix.linearIndex(i,j)] = t*precond(matrix,i,j);
            }
        }
    }

    for (int i = grid.sizeI(); i >= 1; i--)
    {
        for (int j = grid.sizeJ(); j >= 1; j--)
        {
            if (grid.getMaterial(i,j) == FluidCellMaterial::FLUID)
            {
                t = temp[matrix.linearIndex(i,j)] - (matrix.Ax(i,j) * precond(matrix,i,j) * temp[matrix.linearIndex(i+1,j)])
                                                    - (matrix.Ax(i,j) * precond(matrix,i,j) * temp[matrix.linearIndex(i,j+1)]);
                vector[matrix.linearIndex(i,j)] = t*precond(matrix,i,j);
            }
        }
    }
}

void PCGSolver::calcPrecond(const SparseMatrix &matrix, MACFluidGrid &grid)
{
    m_precondCache = DynamicSparseMatrix(1);
    m_precondCache.setGridSize(matrix);
    for (int i = 1; i < grid.sizeI(); i++)
    {
        for (int j = 1; j < grid.sizeJ(); j++)
        {

            if (grid.getMaterial(i,j) == FluidCellMaterial::FLUID)
            {
                precond(matrix,i,j);
            }

        }
    }
}

double PCGSolver::precond(const SparseMatrix &m, int i, int j)
{
    double temp;
    if (m_precondCache.getValue(0,m_precondCache.linearIndex(i,j),temp))
    {
        return temp;
    }

    const double tuning = 0.97;
    const double safety = 0.25;

    temp = m.Adiag(i,j)
            - (m.Ax(i-1,j)*precond(m,i-1,j))*(m.Ax(i-1,j)*precond(m,i-1,j))
            - (m.Ay(i,j-1)*precond(m,i,j-1))*(m.Ay(i,j-1)*precond(m,i,j-1))

            - tuning*(m.Ax(i-1,j)*(m.Ay(i-1,j))*precond(m,i-1,j)*precond(m,i-1,j)
                      + m.Ay(i,j-1)*(m.Ax(i,j-1))*precond(m,i,j-1)*precond(m,i,j-1));

    if(temp < safety*m.Adiag(i,j))
    {
        temp = m.Adiag(i,j);
    }

    temp = 1/std::sqrt(temp);

    m_precondCache.setValue(0,m_precondCache.linearIndex(i,j),temp);

    return temp;
}
