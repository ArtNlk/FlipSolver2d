#ifndef SPARSETRIMATRIX_H
#define SPARSETRIMATRIX_H

#include <vector>
#include <utility>
#include <sstream>

#include "dynamicsparsematrix.h"
#include "linearindexable2d.h"
#include "dynamicuppertriangularsparsematrix.h"

class SparseMatrix : public LinearIndexable2d
{
public:
    typedef std::pair<int,double> StaticRowUnit;

    SparseMatrix(const DynamicSparseMatrix &dynamicMatrix);
    SparseMatrix(const DynamicUpperTriangularSparseMatrix &dynamicMatrix);

    virtual double getValue(int row, int col) const;

    inline int rowCount() const
    {
        return m_rowStart.size()-1;
    };

    inline void getGridSize(int& out_gridSizeI, int& out_gridSizeJ) const
    {
        out_gridSizeI = m_sizeI;
        out_gridSizeJ = m_sizeJ;
    }

    inline double Adiag(int i, int j, MACFluidGrid &grid) const
    {
        ASSERT_BETWEEN(i,-2,m_sizeI);
        ASSERT_BETWEEN(j,-2,m_sizeJ);
        if(i < 0 || j < 0)
        {
            return 0.0;
        }
        int index = grid.linearIndex(i,j);
        return getValue(index,index);
    }

    inline double Ax(int i, int j, MACFluidGrid &grid) const
    {
        ASSERT_BETWEEN(i,-2,m_sizeI);
        ASSERT_BETWEEN(j,-2,m_sizeJ);
        if(i < 0 || j < 0)
        {
            return 0.0;
        }
        int rowIndex = grid.linearIndex(i,j);
        int colIndex = grid.linearIndex(i+1,j);

        return getValue(rowIndex,colIndex);
    }

    inline double Ay(int i, int j, MACFluidGrid &grid) const
    {
        ASSERT_BETWEEN(i,-2,m_sizeI);
        ASSERT_BETWEEN(j,-2,m_sizeJ);
        if(i < 0 || j < 0)
        {
            return 0.0;
        }
        int rowIndex = grid.linearIndex(i,j);
        int colIndex = grid.linearIndex(i,j+1);

        return getValue(rowIndex,colIndex);
    }

    std::vector<double> operator*(const std::vector<double> &v) const;

    inline std::string toString()
    {
        std::ostringstream output;
        for(int i = 0; i < m_sizeI*m_sizeJ; i++)
        {
            output << "|";
            for(int j = 0; j < m_sizeI*m_sizeJ; j++)
            {
                output << "\t" << getValue(i,j) << ",";
            }
            output << "|\n";
        }

        return output.str();
    }

protected:
    std::vector<StaticRowUnit> m_values;
    std::vector<int> m_rowStart;
    int m_size;
};

#endif // SPARSETRIMATRIX_H
