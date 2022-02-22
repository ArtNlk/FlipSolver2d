#ifndef DYNAMICSPARSEMATRIX_H
#define DYNAMICSPARSEMATRIX_H

#include <vector>
#include <utility>

#include "fluidgrid.h"
#include "linearindexable2d.h"

class SparseMatrix;

class DynamicSparseMatrix : public LinearIndexable2d
{
public:
    typedef std::pair<int,double> SparseRowUnit;
    typedef std::vector<SparseRowUnit> SparseRow;

    DynamicSparseMatrix(int size, int avgRowLength = 7);

    DynamicSparseMatrix(MACFluidGrid &grid);

    void resize(int newSize);
    int size() const;

    inline void setGridSize(const MACFluidGrid &grid)
    {
        grid.getSize(m_sizeI, m_sizeJ);
    }

    inline void setGridSize(const DynamicSparseMatrix &matrix)
    {
        matrix.getGridSize(m_sizeI, m_sizeJ);
    }

    void setGridSize(const SparseMatrix &matrix);

    inline void getGridSize(int& out_gridSizeI, int& out_gridSizeJ) const
    {
        out_gridSizeI = m_sizeI;
        out_gridSizeJ = m_sizeJ;
    }

    inline void setAdiag(int i, int j, double value)
    {
        int index = linearIndex(i,j);
        setValue(index, index, value);
    }

    inline void setAx(int i, int j, double value)
    {
        int rowIndex = linearIndex(i,j);
        int colIndex = linearIndex(i+1,j);

        return setValue(rowIndex,colIndex, value);
    }

    inline void setAy(int i, int j, double value)
    {
        int rowIndex = linearIndex(i,j);
        int colIndex = linearIndex(i,j+1);

        return setValue(rowIndex,colIndex, value);
    }

    inline void modifyAdiag(int i, int j, double value)
    {
        int index = linearIndex(i,j);
        setValue(index, index, getValue(index,index) + value);
    }

    inline void modifyAx(int i, int j, double value)
    {
        int rowIndex = linearIndex(i,j);
        int colIndex = linearIndex(i+1,j);

        return setValue(rowIndex,colIndex, getValue(rowIndex, colIndex) + value);
    }

    inline void modifyAy(int i, int j, double value)
    {
        int rowIndex = linearIndex(i,j);
        int colIndex = linearIndex(i,j+1);

        return setValue(rowIndex,colIndex, getValue(rowIndex, colIndex) + value);
    }

    inline int rowSize(int rowIndex) const { return m_rows[rowIndex].size();}

    inline int elementCount() const { return m_elementCount;}

    inline const std::vector<SparseRow> *data() const { return &m_rows;}

    void setValue(int rowIndex, int columnIndex, double value);
    double getValue(int rowIndex, int columnIndex) const;
    bool getValue(int rowIndex, int columnIndex, double &out) const;

protected:
    std::vector<SparseRow> m_rows;
    int m_size;
    int m_elementCount;
};

#endif // DYNAMICSPARSEMATRIX_H
