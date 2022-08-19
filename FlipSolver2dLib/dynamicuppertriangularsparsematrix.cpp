#include "dynamicuppertriangularsparsematrix.h"

#include "mathfuncs.h"
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
    DynamicUpperTriangularSparseMatrix output(grid.cellCount(),7);

    double scale = SimSettings::stepDt() / (SimSettings::density() * SimSettings::dx() * SimSettings::dx());

    for(int i = 0; i < grid.sizeI(); i++)
    {
        for(int j = 0; j < grid.sizeJ(); j++)
        {
            if(grid.isFluid(i,j))
            {
                //X Neighbors
                if(grid.isFluid(i-1,j))
                {
                    output.addToAdiag(i,j,scale, grid);
                }else if(grid.isEmpty(i-1,j))
                {
                    output.addToAdiag(i,j,scale, grid);
                }

                if(grid.isFluid(i+1,j))
                {
                    output.addToAdiag(i,j,scale, grid);
                    output.addToAx(i,j,-scale, grid);
                } else if(grid.isEmpty(i+1,j))
                {
                    output.addToAdiag(i,j,scale, grid);
                }

                //Y Neighbors
                if(grid.isFluid(i,j-1))
                {
                    output.addToAdiag(i,j,scale, grid);
                }else if(grid.isEmpty(i,j-1))
                {
                    output.addToAdiag(i,j,scale, grid);
                }

                if(grid.isFluid(i,j+1))
                {
                    output.addToAdiag(i,j,scale, grid);
                    output.addToAy(i,j,-scale, grid);
                } else if(grid.isEmpty(i,j+1))
                {
                    output.addToAdiag(i,j,scale, grid);
                }
            }
        }
    }
    return output;
}

DynamicUpperTriangularSparseMatrix DynamicUpperTriangularSparseMatrix::forViscosity(MACFluidGrid &grid)
{
    DynamicUpperTriangularSparseMatrix output(grid.cellCount()*2,7);

    float scaleTwoDt = 2*SimSettings::stepDt() / (SimSettings::dx() * SimSettings::dx());
    float scaleTwoDx = SimSettings::stepDt() / (2 * SimSettings::dx() * SimSettings::dx());

    for(int i = 0; i < grid.sizeI(); i++)
    {
        for(int j = 0; j < grid.sizeJ(); j++)
        {
            float fi = static_cast<float>(i);
            float fj = static_cast<float>(j);
            int currLinearIdx = grid.linearIndex(i,j);
            int vBaseIndex = grid.cellCount();
            int uIpOneLinearIdx = grid.linearIndex(i+1,j);
            int vJpOneLinearIdx = vBaseIndex + grid.linearIndex(i,j+1);
            bool uIpOneIsValid = grid.inBounds(i+1,j) && !grid.uVelocityInside(i+1,j);
            bool vJpOneIsValid = grid.inBounds(i,j+1) && !grid.vVelocityInside(i,j+1);
            ///swap uip/vjp indexes with current index
            //U component
            if(uIpOneIsValid)
            {
                output.addTo(currLinearIdx,uIpOneLinearIdx,SimSettings::density());

                if(grid.uVelocityInside(i,j))
                {
                    output.addTo(currLinearIdx,currLinearIdx,-scaleTwoDt * grid.viscosity(i,j));
                    output.addTo(currLinearIdx,uIpOneLinearIdx,scaleTwoDt * grid.viscosity(i,j));
                }

                if(grid.uVelocityInside(i+2,j))
                {
                    output.addTo(currLinearIdx,grid.linearIndex(i+2,j),-scaleTwoDt * grid.viscosity(i+1,j));
                    output.addTo(currLinearIdx,uIpOneLinearIdx,scaleTwoDt * grid.viscosity(i+1,j));
                }

                if(grid.uVelocityInside(i+1,j-1)
                        && grid.vVelocityInside(i+1,j)
                        && grid.vVelocityInside(i,j))
                {
                    float lerpedViscosity = math::lerpCenteredGrid(fi+0.5f,fj-0.5f,grid.viscosityGrid());
                    output.addTo(currLinearIdx,grid.linearIndex(i+1,j-1),-scaleTwoDx * lerpedViscosity);
                    output.addTo(currLinearIdx,vBaseIndex +
                                 grid.linearIndex(i+1,j),scaleTwoDx * lerpedViscosity);
                    output.addTo(currLinearIdx,vBaseIndex +
                                 grid.linearIndex(i,j),-scaleTwoDx * lerpedViscosity);
                    output.addTo(currLinearIdx,uIpOneLinearIdx,scaleTwoDx * lerpedViscosity);
                }

                if(grid.uVelocityInside(i+1,j+1)
                        && grid.vVelocityInside(i+1,j+1)
                        && grid.vVelocityInside(i,j+1))
                {
                    float lerpedViscosity = math::lerpCenteredGrid(fi+0.5f,fj+0.5f,grid.viscosityGrid());
                    output.addTo(currLinearIdx,grid.linearIndex(i+1,j+1),-scaleTwoDx * lerpedViscosity);
                    output.addTo(currLinearIdx,vBaseIndex +
                                 grid.linearIndex(i+1,j+1),-scaleTwoDx * lerpedViscosity);
                    output.addTo(currLinearIdx,vBaseIndex +
                                 grid.linearIndex(i,j+1),scaleTwoDx * lerpedViscosity);
                    output.addTo(currLinearIdx,uIpOneLinearIdx,scaleTwoDx * lerpedViscosity);
                }
            }

            //V component
            if(vJpOneIsValid)
            {
                output.addTo(currLinearIdx,vJpOneLinearIdx,SimSettings::density());

                if(grid.vVelocityInside(i,j))
                {
                    output.addTo(vBaseIndex + currLinearIdx,vBaseIndex + currLinearIdx,
                                 -scaleTwoDt * grid.viscosity(i,j));
                    output.addTo(vBaseIndex + currLinearIdx,vJpOneLinearIdx,scaleTwoDt * grid.viscosity(i,j));
                }

                if(grid.vVelocityInside(i,j+2))
                {
                    output.addTo(vBaseIndex + currLinearIdx,vBaseIndex + grid.linearIndex(i+2,j),
                                 -scaleTwoDt * grid.viscosity(i,j+1));
                    output.addTo(vBaseIndex + currLinearIdx,vJpOneLinearIdx,scaleTwoDt * grid.viscosity(i,j+1));
                }

                if(grid.uVelocityInside(i,j+1)
                        && grid.uVelocityInside(i,j)
                        && grid.vVelocityInside(i-1,j+1))
                {
                    float lerpedViscosity = math::lerpCenteredGrid(fi-0.5f,fj+0.5f,grid.viscosityGrid());
                    output.addTo(vBaseIndex + currLinearIdx,grid.linearIndex(i,j+1),
                                 scaleTwoDx * lerpedViscosity);
                    output.addTo(vBaseIndex + currLinearIdx, grid.linearIndex(i,j),
                                 -scaleTwoDx * lerpedViscosity);
                    output.addTo(vBaseIndex + currLinearIdx,vBaseIndex + grid.linearIndex(i-1,j+1),
                                 -scaleTwoDx * lerpedViscosity);
                    output.addTo(vBaseIndex + currLinearIdx,vJpOneLinearIdx,scaleTwoDx * lerpedViscosity);
                }

                if(grid.uVelocityInside(i+1,j)
                        && !grid.uVelocityInside(i+1,j+1)
                        && !grid.vVelocityInside(i+1,j+1))
                {
                    float lerpedViscosity = math::lerpCenteredGrid(fi+0.5f,fj+0.5f,grid.viscosityGrid());
                    output.addTo(vBaseIndex + currLinearIdx,grid.linearIndex(i+1,j+1),
                                 -scaleTwoDx * lerpedViscosity);
                    output.addTo(vBaseIndex + currLinearIdx, grid.linearIndex(i+1,j),
                                 scaleTwoDx * lerpedViscosity);
                    output.addTo(vBaseIndex + currLinearIdx,vBaseIndex + grid.linearIndex(i+1,j+1),
                                 -scaleTwoDx * lerpedViscosity);
                    output.addTo(vBaseIndex + currLinearIdx,vJpOneLinearIdx,scaleTwoDx * lerpedViscosity);
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
