#include "dynamicuppertriangularsparsematrix.h"

#include "uppertriangularmatrix.h"
#include "simsettings.h"

DynamicUpperTriangularSparseMatrix::DynamicUpperTriangularSparseMatrix(int size, int avgRowLength) :
    LinearIndexable2d(size,size),
    m_rows(size),
    m_size(size),
    m_elementCount(0)
{
    for(int i = 0; i < size; i++)
    {
        m_rows[i].reserve(avgRowLength);
    }
}

DynamicUpperTriangularSparseMatrix DynamicUpperTriangularSparseMatrix::forPressureProjection(MACFluidGrid &grid)
{
    DynamicUpperTriangularSparseMatrix output(grid.fluidCellCount(),7);
    output.m_elementCount = 0;
    for(int i = 0; i < output.m_size; i++)
    {
        output.m_rows[i].reserve(7);
    }
    double scale = SimSettings::stepDt() / (SimSettings::density() * SimSettings::dx() * SimSettings::dx());

    for(int i = 0; i < grid.sizeI(); i++)
    {
        for(int j = 0; j < grid.sizeI(); j++)
        {
            if(grid.isFluid(i,j))
            {
                //X Neighbors
                if(grid.isFluid(i-1,j))
                {
                    output.modifyAdiag(i,j,scale, grid);
                }else if(grid.isEmpty(i-1,j))
                {
                    output.modifyAdiag(i,j,scale, grid);
                }

                if(grid.isFluid(i+1,j))
                {
                    output.modifyAdiag(i,j,scale, grid);
                    output.modifyAx(i,j,-scale, grid);
                } else if(grid.isEmpty(i+1,j))
                {
                    output.modifyAdiag(i,j,scale, grid);
                }

                //Y Neighbors
                if(grid.isFluid(i,j-1))
                {
                    output.modifyAdiag(i,j,scale, grid);
                }else if(grid.isEmpty(i,j-1))
                {
                    output.modifyAdiag(i,j,scale, grid);
                }

                if(grid.isFluid(i,j+1))
                {
                    output.modifyAdiag(i,j,scale, grid);
                    output.modifyAy(i,j,-scale, grid);
                } else if(grid.isEmpty(i,j+1))
                {
                    output.modifyAdiag(i,j,scale, grid);
                }
            }
        }
    }
    return output;
}

void DynamicUpperTriangularSparseMatrix::resize(int newSize)
{
    m_size = newSize;
    m_rows.resize(newSize);
}

int DynamicUpperTriangularSparseMatrix::size() const
{
    return m_size;
}

void DynamicUpperTriangularSparseMatrix::setGridSize(const UpperTriangularMatrix &matrix)
{
    matrix.getGridSize(m_sizeI, m_sizeJ);
}

void DynamicUpperTriangularSparseMatrix::setValue(int rowIndex, int columnIndex, double value)
{
    internalSetValue(rowIndex, columnIndex, value);
    internalSetValue(columnIndex, rowIndex, value);
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

bool DynamicUpperTriangularSparseMatrix::getValue(int rowIndex, int columnIndex, double &out) const
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

    out = 0;
    return false;
}

void DynamicUpperTriangularSparseMatrix::internalSetValue(int rowIndex, int columnIndex, double value)
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
