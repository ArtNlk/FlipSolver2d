#ifndef PCGSOLVER_H
#define PCGSOLVER_H

#include "sparsematrix.h"
#include "dynamicsparsematrix.h"
#include "fluidgrid.h"

class PCGSolver
{
public:
    PCGSolver();

    bool solve(const SparseMatrix &matrix, MACFluidGrid &grid, std::vector<double> &result, const std::vector<double> &vec, int iterLimit = 200);

protected:
    void applyICPrecond(const SparseMatrix &matrix, std::vector<double> &vector, MACFluidGrid &grid);
    void calcPrecond(const SparseMatrix &matrix, MACFluidGrid &grid);
    double precond(const SparseMatrix &matrix, int i, int j);

    std::vector<double> m_residual;
    std::vector<double> m_aux;
    std::vector<double> m_search;
    DynamicSparseMatrix m_precondCache;
    static const double m_tol;
};

#endif // PCGSOLVER_H
