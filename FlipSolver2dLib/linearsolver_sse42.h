#ifndef LINEARSOLVER_SSE42_H
#define LINEARSOLVER_SSE42_H

#include <vector>

#include "materialgrid.h"
#include "threadpool.h"

class LinearSolver;

class LinearSolver_sse42
{
public:
    LinearSolver_sse42() = default;

    static  void dampedJacobiThread(LinearSolver*, const Range range, const MaterialGrid& materials, std::vector<double> &vout,
                                    const std::vector<double> &pressures, const std::vector<double> &rhs);

    static std::pair<std::array<int,5>,std::array<double,5>> getMultigridMatrixEntriesForCell(const MaterialGrid &materials, int linearIdx);
};

#endif // LINEARSOLVER_SSE42_H
