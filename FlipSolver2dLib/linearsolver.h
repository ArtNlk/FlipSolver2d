#ifndef PCGSOLVER_H
#define PCGSOLVER_H

#include <functional>
#include <unordered_map>
#include <vector>

#include "linearindexable2d.h"
#include "materialgrid.h"
#include "threadpool.h"
#include "staticmatrix.h"

#include "pressuredata.h"
#include "PressureIPPCoeficients.h"

class LinearSolver
{
public:
    LinearSolver();
    using SparseMatRowElements = std::array<std::pair<int,double>,5>;

    bool solve(const IndexedPressureParameters &matrixIn,
               const IndexedIPPCoefficients &precond,
               std::vector<double> &result,
               const std::vector<double> &vec,
               int iterLimit = 20,
               double tol = 1e-6);

    friend class LinearSolver_sse42;
};
#endif // PCGSOLVER_H
