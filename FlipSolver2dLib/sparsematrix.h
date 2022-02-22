#ifndef SPARSETRIMATRIX_H
#define SPARSETRIMATRIX_H

#include <vector>
#include <utility>

#include "dynamicsparsematrix.h"
#include "linearindexable2d.h"

class SparseMatrix : public LinearIndexable2d
{
public:
    typedef std::pair<int,double> StaticRowUnit;

    SparseMatrix(const DynamicSparseMatrix &dynamicMatrix);

    double getValue(int row, int col) const;

    inline int rowCount() const
    {
        return m_rowStart.size()-1;
    };

    inline void getGridSize(int& out_gridSizeI, int& out_gridSizeJ) const
    {
        out_gridSizeI = m_sizeI;
        out_gridSizeJ = m_sizeJ;
    }

    inline double Adiag(int i, int j) const
    {
        int index = linearIndex(i,j);
        return getValue(index,index);
    }

    inline double Ax(int i, int j) const
    {
        int rowIndex = linearIndex(i,j);
        int colIndex = linearIndex(i+1,j);

        return getValue(rowIndex,colIndex);
    }

    inline double Ay(int i, int j) const
    {
        int rowIndex = linearIndex(i,j);
        int colIndex = linearIndex(i,j+1);

        return getValue(rowIndex,colIndex);
    }

    std::vector<double> operator*(const std::vector<double> &v) const;

protected:
    std::vector<StaticRowUnit> m_values;
    std::vector<int> m_rowStart;
    int m_size;
};

#endif // SPARSETRIMATRIX_H
