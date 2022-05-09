#include "uppertriangularmatrix.h"

UpperTriangularMatrix::UpperTriangularMatrix(const DynamicUpperTriangularSparseMatrix &dynamicMatrix) :
    SparseMatrix(dynamicMatrix)
{
}

double UpperTriangularMatrix::getValue(int row, int col) const
{
    if(row == -1 || col == -1) {return 0.0;}
    //if(row > col) {std::swap(row,col);}
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
