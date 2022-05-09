#include "sparsematrix.h"

SparseMatrix::SparseMatrix(const DynamicSparseMatrix &dynamicMatrix) :
    LinearIndexable2d(dynamicMatrix.sizeI(), dynamicMatrix.sizeJ()),
    m_rowStart(dynamicMatrix.size()+1),
    m_size(dynamicMatrix.size())
{
    m_values.resize(dynamicMatrix.elementCount());

    m_rowStart[0] = 0;
    m_rowStart[m_size] = dynamicMatrix.elementCount();
    for(int i = 1; i < m_size; i++)
    {
        m_rowStart[i] = m_rowStart[i-1] + dynamicMatrix.rowSize(i-1);
    }

    int totalIndex = 0;
    const std::vector<DynamicSparseMatrix::SparseRow> *rows = dynamicMatrix.data();
    for(int i = 0; i < m_size; i++)
    {
        for(int j = 0; j < (*rows)[i].size(); j++)
        {
            m_values[totalIndex] = (*rows)[i][j];
            totalIndex++;
        }
    }
}

SparseMatrix::SparseMatrix(const DynamicUpperTriangularSparseMatrix &dynamicMatrix):
    LinearIndexable2d(dynamicMatrix.sizeI(), dynamicMatrix.sizeJ()),
    m_rowStart(dynamicMatrix.size()+1),
    m_size(dynamicMatrix.size())
{
    m_values.resize(dynamicMatrix.elementCount());

    m_rowStart[0] = 0;
    m_rowStart[m_size] = dynamicMatrix.elementCount();
    for(int i = 1; i < m_size; i++)
    {
        m_rowStart[i] = m_rowStart[i-1] + dynamicMatrix.rowSize(i-1);
    }

    int totalIndex = 0;
    const std::vector<DynamicSparseMatrix::SparseRow> *rows = dynamicMatrix.data();
    for(int i = 0; i < m_size; i++)
    {
        for(int j = 0; j < (*rows)[i].size(); j++)
        {
            m_values[totalIndex] = (*rows)[i][j];
            totalIndex++;
        }
    }
}

double SparseMatrix::getValue(int row, int col) const
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

std::vector<double> SparseMatrix::operator*(const std::vector<double> &v) const
{
    std::vector<double> output(v.size());
    for(int i = 0; i < m_size; i++)
    {
        output[i] = 0;
        for(int j = m_rowStart[i]; j < m_rowStart[i+1]; j++)
        {
            output[i] += m_values[j].second * v[m_values[j].first];
        }
    }

    return output;
}
