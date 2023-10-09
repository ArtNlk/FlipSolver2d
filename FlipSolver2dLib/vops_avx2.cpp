#include "vops_avx2.h"

void VOps_avx::addMulThread(Range range, std::vector<double> &output, const std::vector<double> &vec1, const std::vector<double> &vec2, double value)
{
    const unsigned int regSize = 4;
    const unsigned int leftover = range.size() % regSize;
    const unsigned int wholeElements = range.size() - leftover;
    const unsigned int newEnd = range.start + wholeElements;
    const unsigned int leftoverEnd = newEnd + leftover;

    __m256d _temp, _res, _val;

    _val = _mm256_set1_pd(value);

    for(unsigned int i = range.start; i < newEnd; i+=regSize)
    {
        //output[i] = vec1[i] + vec2[i]*value;
        _temp = _mm256_loadu_pd(&vec2[i]);
        _temp = _mm256_mul_pd(_temp,_val);
        _res = _mm256_loadu_pd(&vec1[i]);
        _res = _mm256_add_pd(_res,_temp);
        _mm256_storeu_pd(&output[i],_res);
    }

    for(unsigned int i = newEnd; i < leftoverEnd; i++)
    {
        output[i] = vec1[i] + vec2[i]*value;
    }
}
