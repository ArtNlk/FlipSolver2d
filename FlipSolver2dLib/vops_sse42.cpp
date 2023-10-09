#include "vops_sse42.h"

#include <nmmintrin.h>

void VOps_sse42::addMulThread(Range range, std::vector<double> &output, const std::vector<double> &vec1, const std::vector<double> &vec2, double value)
{
    const unsigned int regSize = 2 * 2;
    const unsigned int leftover = range.size() % regSize;
    const unsigned int wholeElements = range.size() - leftover;
    const unsigned int newEnd = range.start + wholeElements;
    const unsigned int leftoverEnd = newEnd + leftover;

    __m128d _temp, _res, _val;
    __m128d _temp1, _res1, _val1;

    _val = _mm_set1_pd(value);

    for(unsigned int i = range.start; i < newEnd; i+=regSize)
    {
        //output[i] = vec1[i] + vec2[i]*value;
        _temp = _mm_loadu_pd(&vec2[i]);
        _temp1 = _mm_loadu_pd(&vec2[i+2]);
        _res = _mm_loadu_pd(&vec1[i]);
        _res1 = _mm_loadu_pd(&vec1[i+2]);

        _temp = _mm_mul_pd(_temp,_val);
        _temp1 = _mm_mul_pd(_temp1,_val);

        _res = _mm_add_pd(_res,_temp);
        _res1 = _mm_add_pd(_res1,_temp1);

        _mm_storeu_pd(&output[i],_res);
        _mm_storeu_pd(&output[i+2],_res1);
    }

    for(unsigned int i = newEnd; i < leftoverEnd; i++)
    {
        output[i] = vec1[i] + vec2[i]*value;
    }
}
