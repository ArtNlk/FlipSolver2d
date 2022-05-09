#ifndef VMATH_H
#define VMATH_H

#include <vector>
#include <cmath>
#include <customassert.h>

namespace vmath
{
inline double dot(std::vector<double> &v1, std::vector<double> &v2)
{
    ASSERT(v1.size() == v2.size())
    double result = 0;
    for(int i = 0; i < v1.size(); i++)
    {
        result += v1[i] * v2[i];
    }

    return result;
}

inline double isZero(const std::vector<double> &v1,const double eps = 1.0e-15)
{
    for(int i = 0; i < v1.size(); i++)
    {
        if(v1[i] > eps)
        {
            return false;
        }
    }

    return true;
}

inline void addMul(std::vector<double> &output,const std::vector<double> &vec1, const std::vector<double> &vec2, double value)
{
    ASSERT(output.size() == vec1.size() && output.size() == vec2.size())
    for(int i = 0; i < output.size(); i++)
    {
        output[i] = vec1[i] + vec2[i]*value;
    }
}

inline void subMul(std::vector<double> &output,const std::vector<double> &vec1, const std::vector<double> &vec2, double value)
{
    ASSERT(output.size() == vec1.size() && output.size() == vec2.size())
    for(int i = 0; i < output.size(); i++)
    {
        output[i] = vec1[i] - vec2[i]*value;
    }
}

inline double maxAbs(std::vector<double> &vec)
{
    int max = 0;
    for(int i = 0; i < vec.size(); i++)
    {
        if(std::abs(vec[i]) > std::abs(vec[max]))
        {
            max = i;
        }
    }

    return std::abs(vec[max]);
}
}
#endif // VMATH_H
