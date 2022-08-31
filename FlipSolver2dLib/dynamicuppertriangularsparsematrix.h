#ifndef DYNAMICUPPERTRIANGULARSPARSEMATRIX_H
#define DYNAMICUPPERTRIANGULARSPARSEMATRIX_H

#include <vector>
#include <utility>
#include <sstream>

#include "fluidgrid.h"
#include "linearindexable2d.h"

class UpperTriangularMatrix;
class Logger;

class DynamicUpperTriangularSparseMatrix : public LinearIndexable2d
{
public:
    typedef std::pair<int,double> SparseRowUnit;
    typedef std::vector<SparseRowUnit> SparseRow;
    friend Logger;

    DynamicUpperTriangularSparseMatrix(int size, int avgRowLength = 7);

    static DynamicUpperTriangularSparseMatrix forPressureProjection(MACFluidGrid &grid);
    static DynamicUpperTriangularSparseMatrix forViscosity(MACFluidGrid &grid);

    void resize(int newSize);
    int size() const;

    inline void setGridSize(const MACFluidGrid &grid)
    {
        grid.getSize(m_sizeI, m_sizeJ);
    }

    inline void setGridSize(const DynamicUpperTriangularSparseMatrix &matrix)
    {
        matrix.getGridSize(m_sizeI, m_sizeJ);
    }

    void setGridSize(const UpperTriangularMatrix &matrix);

    inline void getGridSize(int& out_gridSizeI, int& out_gridSizeJ) const
    {
        out_gridSizeI = m_sizeI;
        out_gridSizeJ = m_sizeJ;
    }

    inline void setAdiag(int i, int j, double value, MACFluidGrid &grid)
    {
        ASSERT_BETWEEN(i,-2,m_sizeI);
        ASSERT_BETWEEN(j,-2,m_sizeJ);
        int index = grid.linearIndex(i,j);
        setValue(index, index, value);
    }

    inline void setAx(int i, int j, double value, MACFluidGrid &grid)
    {
        ASSERT_BETWEEN(i,-2,m_sizeI);
        ASSERT_BETWEEN(j,-2,m_sizeJ);
        int rowIndex = grid.linearIndex(i,j);
        int colIndex = grid.linearIndex(i+1,j);

        return setValue(rowIndex,colIndex, value);
    }

    inline void setAy(int i, int j, double value, MACFluidGrid &grid)
    {
        ASSERT_BETWEEN(i,-2,m_sizeI);
        ASSERT_BETWEEN(j,-2,m_sizeJ);
        int rowIndex = grid.linearIndex(i,j);
        int colIndex = grid.linearIndex(i,j+1);

        return setValue(rowIndex,colIndex, value);
    }

    inline void addTo(int i, int j, double value)
    {
        ASSERT_BETWEEN(i,-2,m_sizeI);
        ASSERT_BETWEEN(j,-2,m_sizeJ);
        setValue(i,j, getValue(i,j) + value);
    }

    inline void addToAdiag(int i, int j, double value, MACFluidGrid &grid)
    {
        ASSERT_BETWEEN(i,-2,m_sizeI);
        ASSERT_BETWEEN(j,-2,m_sizeJ);
        int index = grid.linearIndex(i,j);
        setValue(index, index, getValue(index,index) + value);
    }

    inline void addToAx(int i, int j, double value, MACFluidGrid &grid)
    {
        ASSERT_BETWEEN(i,-2,m_sizeI);
        ASSERT_BETWEEN(j,-2,m_sizeJ);
        int rowIndex = grid.linearIndex(i,j);
        int colIndex = grid.linearIndex(i+1,j);

        return setValue(rowIndex,colIndex, getValue(rowIndex, colIndex) + value);
    }

    inline void addToAy(int i, int j, double value, MACFluidGrid &grid)
    {
        ASSERT_BETWEEN(i,-2,m_sizeI);
        ASSERT_BETWEEN(j,-2,m_sizeJ);
        int rowIndex = grid.linearIndex(i,j);
        int colIndex = grid.linearIndex(i,j+1);

        return setValue(rowIndex,colIndex, getValue(rowIndex, colIndex) + value);
    }

    inline int rowSize(int rowIndex) const { return m_rows[rowIndex].size();}

    inline int elementCount() const { return m_elementCount;}

    inline const std::vector<SparseRow> data() const { return m_rows;}

    void setValue(int rowIndex, int columnIndex, double value);
    double getValue(int rowIndex, int columnIndex) const;
    bool getValue(int rowIndex, int columnIndex, double &out) const;

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

    void internalSetValue(int rowIndex, int columnIndex, double value);

    std::vector<SparseRow> m_rows;
    int m_size;
    int m_elementCount;
};

#endif // DYNAMICUPPERTRIANGULARSPARSEMATRIX_H
