#ifndef STATICMATRIX_H
#define STATICMATRIX_H

#include "dynamicmatrix.h"
#include "linearindexable2d.h"

#include "threadpool.h"

using StaticRowUnit = std::pair<int, double>;

class StaticMatrix
{
public:
    
    StaticMatrix(const DynamicMatrix &dynamicMatrix);

    double getValue(int row, int col) const;

    double Adiag(int i, int j, LinearIndexable2d &indexer) const;

    double Ax(int i, int j, LinearIndexable2d &indexer) const;

    double Ay(int i, int j, LinearIndexable2d &indexer) const;

    std::vector<double> operator*(std::vector<double> &v) const;

    std::string toString();

    size_t size() const;

    void mulThread(Range range, std::vector<double> &vin, std::vector<double>& vout) const;
protected:
    std::vector<StaticRowUnit> m_values;
    std::vector<int> m_rowStart;

    size_t m_size;
};

#endif // STATICMATRIX_H
