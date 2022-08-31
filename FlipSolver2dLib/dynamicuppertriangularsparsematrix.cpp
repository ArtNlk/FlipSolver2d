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
    int validUSamples = grid.validVelocitySampleCountU();
    int validVSamples = grid.validVelocitySampleCountV();
    DynamicUpperTriangularSparseMatrix output(validUSamples + validVSamples,7);

    float scaleTwoDt = 2*SimSettings::stepDt() / (SimSettings::dx() * SimSettings::dx());
    float scaleTwoDx = SimSettings::stepDt() / (2 * SimSettings::dx() * SimSettings::dx());

    for(int i = 0; i < grid.sizeI(); i++)
    {
        for(int j = 0; j < grid.sizeJ(); j++)
        {
            int currLinearIdx = grid.linearViscosityVelocitySampleIndexU(i,j);
            float diag = 0;
            float offdiag = 0;
            if(currLinearIdx != 0)
            {
                float fi = static_cast<float>(i);
                float fj = static_cast<float>(j);
                int vBaseIndex = validVSamples;
                ///swap uip/vjp indexes with current index
                //U component
                output.addTo(currLinearIdx,currLinearIdx,SimSettings::density());
                diag += std::abs(SimSettings::density());

                int uImOneLinearIdx = grid.linearViscosityVelocitySampleIndexU(i-1,j);

                if(uImOneLinearIdx != -1)
                {
                    output.addTo(currLinearIdx,
                                 uImOneLinearIdx,
                                 -scaleTwoDt * grid.viscosity(i-1,j));

                    output.addTo(currLinearIdx,
                                 currLinearIdx,
                                 scaleTwoDt * grid.viscosity(i-1,j));

                    diag += scaleTwoDt * grid.viscosity(i-1,j);
                    offdiag += std::abs(-scaleTwoDt * grid.viscosity(i-1,j));
                }

                int uIpOneLinearIdx = grid.linearViscosityVelocitySampleIndexU(i+1,j);

                if(uIpOneLinearIdx != 1)
                {
                    output.addTo(currLinearIdx,
                                 grid.linearViscosityVelocitySampleIndexU(i+1,j),
                                 -scaleTwoDt * grid.viscosity(i,j));

                    output.addTo(currLinearIdx,
                                 currLinearIdx,
                                 scaleTwoDt * grid.viscosity(i,j));
                    diag += scaleTwoDt * grid.viscosity(i,j);
                    offdiag += std::abs(-scaleTwoDt * grid.viscosity(i,j));
                }

                if(grid.uVelocityInside(i,j-1)
                        && grid.vVelocityInside(i,j)
                        && grid.vVelocityInside(i-1,j))
                {
                    float lerpedViscosity = math::lerpCenteredGrid(fi-0.5f,fj-0.5f,grid.viscosityGrid());
                    output.addTo(currLinearIdx,
                                 grid.linearViscosityVelocitySampleIndexU(i,j-1),
                                 -scaleTwoDx * lerpedViscosity);

                    output.addTo(currLinearIdx,
                                 vBaseIndex + grid.linearViscosityVelocitySampleIndexV(i,j),
                                 scaleTwoDx * lerpedViscosity);

                    output.addTo(currLinearIdx,
                                 vBaseIndex + grid.linearViscosityVelocitySampleIndexV(i - 1,j),
                                 -scaleTwoDx * lerpedViscosity);

                    output.addTo(currLinearIdx,currLinearIdx,scaleTwoDx * lerpedViscosity);

                    diag += scaleTwoDx * lerpedViscosity;
                    offdiag += std::abs(-scaleTwoDx * lerpedViscosity) * 3;
                }

                if(grid.uVelocityInside(i,j+1)
                        && grid.vVelocityInside(i,j+1)
                        && grid.vVelocityInside(i-1,j+1))
                {
                    float lerpedViscosity = math::lerpCenteredGrid(fi-0.5f,fj+0.5f,grid.viscosityGrid());
                    output.addTo(currLinearIdx,
                                 grid.linearViscosityVelocitySampleIndexU(i,j+1),
                                 -scaleTwoDx * lerpedViscosity);

                    output.addTo(currLinearIdx,
                                 vBaseIndex + grid.linearViscosityVelocitySampleIndexV(i,j+1),
                                 -scaleTwoDx * lerpedViscosity);

                    output.addTo(currLinearIdx,
                                 vBaseIndex + grid.linearViscosityVelocitySampleIndexV(i-1,j+1),
                                 scaleTwoDx * lerpedViscosity);

                    output.addTo(currLinearIdx,currLinearIdx,scaleTwoDx * lerpedViscosity);

                    diag += scaleTwoDx * lerpedViscosity;
                    offdiag += std::abs(scaleTwoDx * lerpedViscosity) * 3;
                }

                if (std::abs(diag) <= offdiag)
                {
                    std::cout << "Non dominant row for U: " << i << ',' << j
                              << "d/offd: " << std::abs(diag) << ',' << offdiag << '\n';
                }
            }

            diag = 0;
            offdiag = 0;

            currLinearIdx = grid.linearViscosityVelocitySampleIndexV(i,j);
            //V component
            if(currLinearIdx != 0)
            {
                float fi = static_cast<float>(i);
                float fj = static_cast<float>(j);
                int vBaseIndex = validVSamples;
                bool currentUIsValid = grid.inBounds(i,j) && grid.uVelocityInside(i,j);
                bool currentVIsValid = grid.inBounds(i,j) && grid.vVelocityInside(i,j);
                output.addTo(vBaseIndex + currLinearIdx,vBaseIndex + currLinearIdx,SimSettings::density());
                diag += SimSettings::density();

                if(grid.vVelocityInside(i,j-1))
                {
                    output.addTo(vBaseIndex + currLinearIdx,
                                 vBaseIndex + grid.linearViscosityVelocitySampleIndexV(i,j-1),
                                 -scaleTwoDt * grid.viscosity(i,j-1));

                    output.addTo(vBaseIndex + currLinearIdx,
                                 vBaseIndex + currLinearIdx,
                                 scaleTwoDt * grid.viscosity(i,j-1));
                    diag += scaleTwoDt * grid.viscosity(i,j-1);
                    offdiag += std::abs(scaleTwoDt * grid.viscosity(i,j-1));
                }

                if(grid.vVelocityInside(i,j+1))
                {
                    output.addTo(vBaseIndex + currLinearIdx,
                                 vBaseIndex + grid.linearViscosityVelocitySampleIndexV(i,j+1),
                                 -scaleTwoDt * grid.viscosity(i,j));

                    output.addTo(vBaseIndex + currLinearIdx,
                                 vBaseIndex + currLinearIdx,
                                 scaleTwoDt * grid.viscosity(i,j));

                    diag += scaleTwoDt * grid.viscosity(i,j);
                    offdiag += std::abs(scaleTwoDt * grid.viscosity(i,j));
                }

                if(grid.uVelocityInside(i,j)
                        && grid.uVelocityInside(i,j-1)
                        && grid.vVelocityInside(i-1,j))
                {
                    float lerpedViscosity = math::lerpCenteredGrid(fi-0.5f,fj-0.5f,grid.viscosityGrid());

                    output.addTo(vBaseIndex + currLinearIdx,
                                 grid.linearViscosityVelocitySampleIndexU(i,j),
                                 scaleTwoDx * lerpedViscosity);

                    output.addTo(vBaseIndex + currLinearIdx,
                                 grid.linearViscosityVelocitySampleIndexU(i,j-1),
                                 -scaleTwoDx * lerpedViscosity);

                    output.addTo(vBaseIndex + currLinearIdx,
                                 vBaseIndex + grid.linearViscosityVelocitySampleIndexV(i-1,j),
                                 -scaleTwoDx * lerpedViscosity);

                    output.addTo(vBaseIndex + currLinearIdx,
                                 vBaseIndex + currLinearIdx,
                                 scaleTwoDx * lerpedViscosity);

                    diag += scaleTwoDx * lerpedViscosity;
                    offdiag += std::abs(scaleTwoDx * lerpedViscosity)*3;
                }

                if(grid.uVelocityInside(i+1,j)
                        && grid.uVelocityInside(i+1,j-1)
                        && grid.vVelocityInside(i+1,j))
                {
                    float lerpedViscosity = math::lerpCenteredGrid(fi+0.5f,fj-0.5f,grid.viscosityGrid());

                    output.addTo(vBaseIndex + currLinearIdx,
                                 grid.linearViscosityVelocitySampleIndexU(i+1,j),
                                 -scaleTwoDx * lerpedViscosity);

                    output.addTo(vBaseIndex + currLinearIdx,
                                 grid.linearViscosityVelocitySampleIndexU(i+1,j-1),
                                 scaleTwoDx * lerpedViscosity);

                    output.addTo(vBaseIndex + currLinearIdx,
                                 vBaseIndex + grid.linearViscosityVelocitySampleIndexV(i+1,j),
                                 -scaleTwoDx * lerpedViscosity);

                    output.addTo(vBaseIndex + currLinearIdx,
                                 vBaseIndex + currLinearIdx,
                                 scaleTwoDx * lerpedViscosity);

                    diag += scaleTwoDx * lerpedViscosity;
                    offdiag += std::abs(scaleTwoDx * lerpedViscosity)*3;
                }

                if (std::abs(diag) <= offdiag)
                {
                    std::cout << "Non dominant row for V: " << i << ',' << j
                              << "d/offd: " << std::abs(diag) << ',' << offdiag << '\n';
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
