#include "dynamicsparsematrix.h"

#include "sparsematrix.h"

DynamicSparseMatrix::DynamicSparseMatrix(int size, int avgRowLength) :
    LinearIndexable3d(-1,-1,-1),
    m_size(size),
    m_rows(size),
    m_elementCount(0)
{
    for(int i = 0; i < size; i++)
    {
        m_rows[i].reserve(avgRowLength);
    }
}

DynamicSparseMatrix::DynamicSparseMatrix(FluidGrid &grid) :
    LinearIndexable3d(grid.sizeI(),grid.sizeJ(),grid.sizeK()),
    m_size(grid.sizeI() * grid.sizeJ() * grid.sizeK()),
    m_rows(grid.sizeI() * grid.sizeJ() * grid.sizeK()),
    m_elementCount(0)
{
    for(int i = 0; i < m_size; i++)
    {
        m_rows[i].reserve(7);
    }
    double scale = grid.getDt() / (grid.getFluidDensity() * grid.getSideLength() * grid.getSideLength());

    for(int i = 0; i < m_sizeI; i++)
    {
        for(int j = 0; j < m_sizeJ; j++)
        {
            for(int k = 0; k < m_sizeK; k++)
            {
                if(grid.getMaterial(i,j,k) == FluidCell::FLUID)
                {
                    //X Neighbors
                    if(grid.getMaterial(i-1,j,k) == FluidCell::FLUID)
                    {
                        modifyAdiag(i,j,k,scale);
                    }
                    if(grid.getMaterial(i+1,j,k) == FluidCell::FLUID)
                    {
                        modifyAdiag(i,j,k,scale);
                        modifyAx(i,j,k,-scale);
                    } else if(grid.getMaterial(i+1,j,k) == FluidCell::AIR)
                    {
                        modifyAdiag(i,j,k,scale);
                    }

                    //Y Neighbors
                    if(grid.getMaterial(i,j-1,k) == FluidCell::FLUID)
                    {
                        modifyAdiag(i,j,k,scale);
                    }
                    if(grid.getMaterial(i,j+1,k) == FluidCell::FLUID)
                    {
                        modifyAdiag(i,j,k,scale);
                        modifyAy(i,j,k,-scale);
                    } else if(grid.getMaterial(i,j+1,k) == FluidCell::AIR)
                    {
                        modifyAdiag(i,j,k,scale);
                    }

                    //Z Neighbors
                    if(grid.getMaterial(i,j,k-1) == FluidCell::FLUID)
                    {
                        modifyAdiag(i,j,k,scale);
                    }
                    if(grid.getMaterial(i,j,k+1) == FluidCell::FLUID)
                    {
                        modifyAdiag(i,j,k,scale);
                        modifyAz(i,j,k,-scale);
                    } else if(grid.getMaterial(i,j,k+1) == FluidCell::AIR)
                    {
                        modifyAdiag(i,j,k,scale);
                    }
                }
            }
        }
    }
}

void DynamicSparseMatrix::resize(int newSize)
{
    m_size = newSize;
    m_rows.resize(newSize);
}

int DynamicSparseMatrix::size() const
{
    return m_size;
}

void DynamicSparseMatrix::setGridSize(const SparseMatrix &matrix)
{
    matrix.getGridSize(m_sizeI, m_sizeJ, m_sizeK);
}

void DynamicSparseMatrix::setValue(int rowIndex, int columnIndex, double value)
{
    SparseRow &targetRow = m_rows[rowIndex];
    for(int column = 0; column < targetRow.size(); column++)
    {
        if(targetRow[column].first == columnIndex)
        {
            targetRow[column].second = value;
            return;
        }
        else if(targetRow[column].first > columnIndex)
        {
            targetRow.insert(targetRow.begin()+column,SparseRowUnit(column,value));
            m_elementCount++;
            return;
        }
    }
    targetRow.push_back(SparseRowUnit(columnIndex,value));
    m_elementCount++;
}

double DynamicSparseMatrix::getValue(int rowIndex, int columnIndex) const
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

bool DynamicSparseMatrix::getValue(int rowIndex, int columnIndex, double &out) const
{
    const SparseRow &targetRow = m_rows[rowIndex];
    for(int column = 0; column < targetRow.size(); column++)
    {
        if(targetRow[column].first == columnIndex)
        {
            out = targetRow[column].second;
            return true;
        }
    }

    return false;
}
