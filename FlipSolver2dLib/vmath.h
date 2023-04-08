#ifndef VMATH_H
#define VMATH_H

#include "threadpool.h"
#include <limits>
#include <vector>
#include <cmath>
#include <customassert.h>

namespace vsimmath
{

void dotThread(Range range, std::vector<double> &v1, std::vector<double> &v2, double& output)
{
    ASSERT(v1.size() == v2.size());
    double result = 0.0;
    for(int i = 0; i < v1.size(); i++)
    {
        result += v1[i] * v2[i];
    }
    output = result;
}

inline double dot(std::vector<double> &v1, std::vector<double> &v2)
{
    ASSERT(v1.size() == v2.size());
    double result = 0;
    std::vector<Range> ranges = ThreadPool::i()->splitRange(v1.size());
    std::vector<double> outputs(ranges.size());

    int i = 0;
    for(Range& range : ranges)
    {
        ThreadPool::i()->enqueue(&dotThread,range,std::ref(v1),std::ref(v2),std::ref(outputs[i]));
        i++;
    }
    ThreadPool::i()->wait();

    for(double& v : outputs)
    {
        result += v;
    }

    return result;
}

inline bool isZero(const std::vector<double> &v1,const double eps = 1.0e-15)
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

void addMulThread(Range range, std::vector<double> &output,const std::vector<double> &vec1, const std::vector<double> &vec2, double value)
{
    for(int i = range.start; i < range.end; i++)
    {
        output[i] = vec1[i] + vec2[i]*value;
    }
}

inline void addMul(std::vector<double> &output,const std::vector<double> &vec1, const std::vector<double> &vec2, double value)
{
    ASSERT(output.size() == vec1.size() && output.size() == vec2.size());
    std::vector<Range> ranges = ThreadPool::i()->splitRange(vec1.size());

    for(Range& range : ranges)
    {
        ThreadPool::i()->enqueue(&addMulThread,range,std::ref(output),std::ref(vec1),std::ref(vec2),value);
    }
    ThreadPool::i()->wait();
}

void subMulThread(Range range, std::vector<double> &output,const std::vector<double> &vec1, const std::vector<double> &vec2, double value)
{
    for(int i = range.start; i < range.end; i++)
    {
        output[i] = vec1[i] - vec2[i]*value;
    }
}

inline void subMul(std::vector<double> &output,const std::vector<double> &vec1, const std::vector<double> &vec2, double value)
{
    ASSERT(output.size() == vec1.size() && output.size() == vec2.size());
    std::vector<Range> ranges = ThreadPool::i()->splitRange(vec1.size());

    for(Range& range : ranges)
    {
        ThreadPool::i()->enqueue(&subMulThread,range,std::ref(output),std::ref(vec1),std::ref(vec2),value);
    }
    ThreadPool::i()->wait();
}

void maxAbsThread(Range range, std::vector<double> &vec, double& max)
{
    double maxAbs = std::numeric_limits<double>::min();
    for(int i = range.start; i < range.end; i++)
    {
        if(std::abs(vec[i]) > max)
        {
            maxAbs = std::abs(vec[i]);
        }
    }
    max = maxAbs;
}

inline double maxAbs(std::vector<double> &vec)
{
    std::vector<Range> ranges = ThreadPool::i()->splitRange(vec.size());
    std::vector<double> maxAbs(ranges.size());
    double max = 0;

    int i = 0;
    for(Range& range : ranges)
    {
        ThreadPool::i()->enqueue(&maxAbsThread,range,std::ref(vec),std::ref(maxAbs[i]));
        i++;
    }
    ThreadPool::i()->wait();

    for(double& v : maxAbs)
    {
        if(v > max)
        {
            max = std::abs(v);
        }
    }

    return max;


//    for(int i = 0; i < vec.size(); i++)
//    {
//        if(std::abs(vec[i]) > std::abs(vec[max]))
//        {
//            max = i;
//        }
//    }

//    return std::abs(vec[max]);
}
}
#endif // VMATH_H
