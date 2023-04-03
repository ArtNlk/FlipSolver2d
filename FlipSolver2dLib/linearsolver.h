#ifndef PCGSOLVER_H
#define PCGSOLVER_H

#include <functional>
#include <unordered_map>

#include "uppertriangularmatrix.h"

class LinearSolver
{
public:
    using SparseMatRowElements = std::array<std::pair<int,double>,5>;
    using MatElementProvider = std::function<SparseMatRowElements(int)>;

    static bool solve(const DynamicUpperTriangularSparseMatrix &matrixIn, std::vector<double> &result, const std::vector<double> &vec, int iterLimit = 20);
    static bool mfcgSolve(MatElementProvider elementProvider, std::vector<double> &result, const std::vector<double> &vec, int iterLimit = 20);
    static void nomatVMul(MatElementProvider elementProvider, const std::vector<double> &vin, std::vector<double> &vout);
protected:
    static void applyICPrecond(const DynamicUpperTriangularSparseMatrix &precond, std::vector<double> const &in, std::vector<double> &out);
    static DynamicUpperTriangularSparseMatrix calcPrecond(const DynamicUpperTriangularSparseMatrix &matrix);
    static const double m_tol;
};

#endif // PCGSOLVER_H
