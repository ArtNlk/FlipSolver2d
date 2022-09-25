#ifndef PCGSOLVER_H
#define PCGSOLVER_H

#include <unordered_map>

#include "uppertriangularmatrix.h"
#include "fluidgrid.h"

class PCGSolver
{
public:
    static bool solve(const DynamicUpperTriangularSparseMatrix &matrixIn, std::vector<double> &result, const std::vector<double> &vec, int iterLimit = 20);

protected:
    static void applyICPrecond(const DynamicUpperTriangularSparseMatrix &precond, std::vector<double> const &in, std::vector<double> &out);
    static DynamicUpperTriangularSparseMatrix calcPrecond(const DynamicUpperTriangularSparseMatrix &matrix);
    static const double m_tol;
};

#endif // PCGSOLVER_H
