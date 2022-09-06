#include "uppertriangularmatrix.h"

#include <stdexcept>

UpperTriangularMatrix::UpperTriangularMatrix(DynamicUpperTriangularSparseMatrix &dynamicMatrix):
    SquareMatrix(dynamicMatrix.size()),
    m_rowStart(dynamicMatrix.size()+1)
{
    m_values.resize(dynamicMatrix.elementCount());

    m_rowStart[0] = 0;
    m_rowStart[size()] = dynamicMatrix.elementCount();
    for(int i = 1; i < size(); i++)
    {
        m_rowStart[i] = m_rowStart[i-1] + dynamicMatrix.rowSize(i-1);
    }

    int totalIndex = 0;
      std::vector<DynamicUpperTriangularSparseMatrix::SparseRow> rows = dynamicMatrix.data();
    for(int i = 0; i < size(); i++)
    {
        for(int j = 0; j < rows[i].size(); j++)
        {
            m_values[totalIndex] = rows[i][j];
            totalIndex++;
        }
    }
}

double UpperTriangularMatrix::getValue(int row, int col) const
{
    if(row == -1 || col == -1) {return 0.0;}
    int rowStartIndex = m_rowStart[row];
    int rowEndIndex = m_rowStart[row+1];
    for(int i = rowStartIndex; i < rowEndIndex; i++)
    {
        if(m_values[i].first == col)
        {
            return m_values[i].second;
        }
    }

    return 0;
}

void UpperTriangularMatrix::setValue(int row, int col, double value)
{
    throw std::runtime_error("Set value is not supported in non-dynamic matrices");
}

int UpperTriangularMatrix::rowSize(int rowIndex)
{
    throw std::runtime_error("Row size is not supported in non-dynamic matrices");
}

double UpperTriangularMatrix::Adiag(int i, int j, MACFluidGrid &grid) const
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

double UpperTriangularMatrix::Ax(int i, int j, MACFluidGrid &grid) const
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

double UpperTriangularMatrix::Ay(int i, int j, MACFluidGrid &grid) const
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

std::vector<double> UpperTriangularMatrix::operator*(std::vector<double> &v) const
{
    std::vector<double> output(v.size());
    for(int i = 0; i < size(); i++)
    {
        output[i] = 0;
        for(int j = m_rowStart[i]; j < m_rowStart[i+1]; j++)
        {
            output[i] += m_values[j].second * v[m_values[j].first];
        }
    }

    return output;
}

std::string UpperTriangularMatrix::toString()
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
