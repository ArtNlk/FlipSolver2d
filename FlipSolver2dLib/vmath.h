#ifndef VMATH_H
#define VMATH_H

#include "threadpool.h"
#include <limits>
#include <vector>
#include <cmath>
#include <customassert.h>

class VOps
{
public:
    static VOps& i();
    VOps(VOps &other) = delete;
    void operator=(const VOps &) = delete;

    double dot(std::vector<double> &v1, std::vector<double> &v2);

    bool isZero(const std::vector<double> &v1,const double eps = 1.0e-15);

    void addMul(std::vector<double> &output,const std::vector<double> &vec1, const std::vector<double> &vec2, double value);

    void subMul(std::vector<double> &output,const std::vector<double> &vec1, const std::vector<double> &vec2, double value);

    double maxAbs(std::vector<double> &vec);

protected:
    VOps();

    static void dotThread(Range range, std::vector<double> &v1, std::vector<double> &v2, double& output);

    static void addMulThread(Range range, std::vector<double> &output,const std::vector<double> &vec1,
                      const std::vector<double> &vec2, double value);

    static void subMulThread(Range range, std::vector<double> &output,const std::vector<double> &vec1,
                      const std::vector<double> &vec2, double value);

    static void maxAbsThread(Range range, std::vector<double> &vec, double& max);

    void (*m_dotThread)(Range,std::vector<double>&, std::vector<double>&, double&);
    void (*m_addMulThread)(Range, std::vector<double>&,const std::vector<double>&,const std::vector<double>&,double);
    void (*m_subMulThread)(Range range, std::vector<double>&,const std::vector<double>&,const std::vector<double>&,double);
    void (*m_maxAbsThread)(Range range, std::vector<double>&,double&);
};
#endif // VMATH_H
