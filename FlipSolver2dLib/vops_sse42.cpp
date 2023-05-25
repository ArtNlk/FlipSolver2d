#include "vops_sse42.h"

#include <nmmintrin.h>

void VOps_sse42::addMulThread(Range range, std::vector<double> &output, const std::vector<double> &vec1, const std::vector<double> &vec2, double value)
{
    const unsigned int regSize = 2;
    const unsigned int leftover = range.size() % regSize;
    const unsigned int wholeElements = range.size() - leftover;
    const unsigned int newEnd = range.start + wholeElements;
    const unsigned int leftoverEnd = newEnd + leftover;

    __m128d _temp, _res, _val;

    _val = _mm_set1_pd(value);

    for(unsigned int i = range.start; i < newEnd; i+=regSize)
    {
        //output[i] = vec1[i] + vec2[i]*value;
        _temp = _mm_set_pd(vec2[i+1],vec2[i]);
        _temp = _mm_mul_pd(_temp,_val);
        _res = _mm_set_pd(vec1[i+1],vec1[i]);
        _res = _mm_add_pd(_res,_temp);
        output[i] = _res[0];
        output[i+1] = _res[1];
    }

    for(unsigned int i = newEnd; i < leftoverEnd; i++)
    {
        output[i] = vec1[i] + vec2[i]*value;
    }
}
