#ifndef PCGSOLVER_H
#define PCGSOLVER_H

#include <vector>

#include "pressuredata.h"
#include "PressureIPPCoeficients.h"

class LinearSolver
{
public:
    LinearSolver() = default;
    using SparseMatRowElements = std::array<std::pair<int,double>,5>;

    int solve(const IndexedPressureParameters &matrixIn,
               const IndexedIPPCoefficients &precond,
               std::vector<double> &result,
               const std::vector<double> &vec,
               int iterLimit = 20,
               double tol = 1e-6);

    friend class LinearSolver_sse42;
};
#endif // PCGSOLVER_H
