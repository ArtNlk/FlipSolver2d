#include "dynamicsparsematrix.h"

#include "sparsematrix.h"
#include "simsettings.h"

DynamicSparseMatrix::DynamicSparseMatrix(int size, int avgRowLength) :
    LinearIndexable2d(-1,-1),
    m_rows(size),
    m_size(size),
    m_elementCount(0)
{
    for(int i = 0; i < size; i++)
    {
        m_rows[i].reserve(avgRowLength);
    }
}

DynamicSparseMatrix::DynamicSparseMatrix(MACFluidGrid &grid) :
    LinearIndexable2d(grid.sizeI(),grid.sizeJ()),
    m_rows(grid.sizeI() * grid.sizeJ()),
    m_size(grid.sizeI() * grid.sizeJ()),
    m_elementCount(0)
{
    for(int i = 0; i < m_size; i++)
    {
        m_rows[i].reserve(7);
    }
    double scale = SimSettings::dt() / (SimSettings::density() * SimSettings::dx() * SimSettings::dx());

    for(int i = 0; i < m_sizeI; i++)
    {
        for(int j = 0; j < m_sizeJ; j++)
        {
            if(grid.getMaterial(i,j) == FluidCellMaterial::FLUID)
            {
                //X Neighbors
                if(grid.getMaterial(i-1,j) == FluidCellMaterial::FLUID)
                {
                    modifyAdiag(i,j,scale);
                }
                if(grid.getMaterial(i+1,j) == FluidCellMaterial::FLUID)
                {
                    modifyAdiag(i,j,scale);
                    modifyAx(i,j,-scale);
                } else if(grid.getMaterial(i+1,j) == FluidCellMaterial::EMPTY)
                {
                    modifyAdiag(i,j,scale);
                }

                //Y Neighbors
                if(grid.getMaterial(i,j-1) == FluidCellMaterial::FLUID)
                {
                    modifyAdiag(i,j,scale);
                }
                if(grid.getMaterial(i,j+1) == FluidCellMaterial::FLUID)
                {
                    modifyAdiag(i,j,scale);
                    modifyAy(i,j,-scale);
                } else if(grid.getMaterial(i,j+1) == FluidCellMaterial::EMPTY)
                {
                    modifyAdiag(i,j,scale);
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
    matrix.getGridSize(m_sizeI, m_sizeJ);
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
