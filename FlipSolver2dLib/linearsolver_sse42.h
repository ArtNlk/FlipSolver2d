#ifndef LINEARSOLVER_SSE42_H
#define LINEARSOLVER_SSE42_H

#include <vector>

#include "materialgrid.h"
#include "threadpool.h"

class LinearSolver_sse42
{
public:
    LinearSolver_sse42() = default;

    static  void premaskPressuresThread(const Range range, const MaterialGrid& materials,
                                       std::vector<double> &pressures);
};

#endif // LINEARSOLVER_SSE42_H
