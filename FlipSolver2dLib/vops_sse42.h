#ifndef VOPS_SSE42_H
#define VOPS_SSE42_H

#include <vector>

#include "threadpool.h"

class VOps_sse42
{
public:
    VOps_sse42() = default;

    static void addMulThread(Range range, std::vector<double> &output,const std::vector<double> &vec1,
                             const std::vector<double> &vec2, double value);
};

#endif // VOPS_SSE42_H
