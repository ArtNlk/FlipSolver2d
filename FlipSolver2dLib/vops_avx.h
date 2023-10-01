#ifndef VOPS_AVX_H
#define VOPS_AVX_H

#include <immintrin.h>

#include "threadpool.h"
class VOps_avx
{
public:
    VOps_avx();

    static void addMulThread(Range range, std::vector<double> &output,const std::vector<double> &vec1,
                             const std::vector<double> &vec2, double value);
};

#endif // VOPS_AVX_H
