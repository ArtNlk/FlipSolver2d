#include "dynamicuppertriangularsparsematrix.h"

#include "mathfuncs.h"
#include "simsettings.h"

DynamicUpperTriangularSparseMatrix::DynamicUpperTriangularSparseMatrix(int size, int avgRowLength) :
    SquareMatrix(size),
    m_rows(size),
    m_size(size),
    m_elementCount(0)
{
    for(int i = 0; i < size; i++)
    {
        m_rows[i].reserve(avgRowLength);
    }
}


void DynamicUpperTriangularSparseMatrix::setAdiag(int i, int j, double value, MACFluidGrid &grid)
{
    ASSERT_BETWEEN(i,-2,m_sizeI);
    ASSERT_BETWEEN(j,-2,m_sizeJ);
    int index = grid.linearIndex(i,j);
    setValue(index, index, value);
}

void DynamicUpperTriangularSparseMatrix::setAx(int i, int j, double value, MACFluidGrid &grid)
{
    ASSERT_BETWEEN(i,-2,m_sizeI);
    ASSERT_BETWEEN(j,-2,m_sizeJ);
    int rowIndex = grid.linearIndex(i,j);
    int colIndex = grid.linearIndex(i+1,j);

    setValue(rowIndex,colIndex, value);
    setValue(colIndex,rowIndex, value);
}

void DynamicUpperTriangularSparseMatrix::setAy(int i, int j, double value, MACFluidGrid &grid)
{
    ASSERT_BETWEEN(i,-2,m_sizeI);
    ASSERT_BETWEEN(j,-2,m_sizeJ);
    int rowIndex = grid.linearIndex(i,j);
    int colIndex = grid.linearIndex(i,j+1);

    setValue(rowIndex,colIndex, value);
    setValue(colIndex,rowIndex, value);
}

void DynamicUpperTriangularSparseMatrix::addTo(int i, int j, double value)
{
    ASSERT_BETWEEN(i,-2,m_sizeI);
    ASSERT_BETWEEN(j,-2,m_sizeJ);
    setValue(i,j, getValue(i,j) + value);
}

void DynamicUpperTriangularSparseMatrix::addToAdiag(int i, int j, double value, MACFluidGrid &grid)
{
    ASSERT_BETWEEN(i,-2,m_sizeI);
    ASSERT_BETWEEN(j,-2,m_sizeJ);
    int index = grid.linearIndex(i,j);
    setValue(index, index, getValue(index,index) + value);
}

void DynamicUpperTriangularSparseMatrix::addToAx(int i, int j, double value, MACFluidGrid &grid)
{
    ASSERT_BETWEEN(i,-2,m_sizeI);
    ASSERT_BETWEEN(j,-2,m_sizeJ);
    int rowIndex = grid.linearIndex(i,j);
    int colIndex = grid.linearIndex(i+1,j);

    return setValue(rowIndex,colIndex, getValue(rowIndex, colIndex) + value);
}

void DynamicUpperTriangularSparseMatrix::addToAy(int i, int j, double value, MACFluidGrid &grid)
{
    ASSERT_BETWEEN(i,-2,m_sizeI);
    ASSERT_BETWEEN(j,-2,m_sizeJ);
    int rowIndex = grid.linearIndex(i,j);
    int colIndex = grid.linearIndex(i,j+1);

    return setValue(rowIndex,colIndex, getValue(rowIndex, colIndex) + value);
}

int DynamicUpperTriangularSparseMatrix::rowSize(int rowIndex) { return m_rows[rowIndex].size();}

int DynamicUpperTriangularSparseMatrix::elementCount() { return m_elementCount;}

const std::vector<DynamicUpperTriangularSparseMatrix::SparseRow> DynamicUpperTriangularSparseMatrix::data() const { return m_rows;}

void DynamicUpperTriangularSparseMatrix::setValue(int rowIndex, int columnIndex, double value)
{
    SparseRow &targetRow = m_rows[rowIndex];
    for(int i = 0; i < targetRow.size(); i++)
    {
        if(targetRow[i].first == columnIndex)
        {
            targetRow[i].second = value;
            return;
        }
        else if(targetRow[i].first > columnIndex)
        {
            targetRow.insert(targetRow.begin()+i,SparseRowUnit(columnIndex,value));
            m_elementCount++;
            return;
        }
    }
    targetRow.push_back(SparseRowUnit(columnIndex,value));
    m_elementCount++;
}

double DynamicUpperTriangularSparseMatrix::getValue(int rowIndex, int columnIndex) const
{
    const SparseRow &targetRow = m_rows[rowIndex];
    for(int column = 0; column < targetRow.size(); column++)
    {
        if(targetRow[column].first == columnIndex)
        {
            return targetRow[column].second;
        }
    }

    return 0;
}

std::string DynamicUpperTriangularSparseMatrix::toString()
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
