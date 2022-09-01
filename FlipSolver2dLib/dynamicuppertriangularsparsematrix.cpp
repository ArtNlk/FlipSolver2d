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
            int currLinearIdxU = grid.linearViscosityVelocitySampleIndexU(i,j);
            int currLinearIdxV = grid.linearViscosityVelocitySampleIndexV(i,j);
            float diag = 0;
            float offdiag = 0;
            if(currLinearIdxU != -1)
            {
                float fi = static_cast<float>(i);
                float fj = static_cast<float>(j);
                int vBaseIndex = validUSamples;
                ///swap uip/vjp indexes with current index
                //U component
                output.addTo(currLinearIdxU,currLinearIdxU,SimSettings::density());
                diag += std::abs(SimSettings::density());

                int uImOneLinearIdx = grid.linearViscosityVelocitySampleIndexU(i-1,j);

                if(uImOneLinearIdx != -1)
                {
                    output.addTo(currLinearIdxU,
                                 uImOneLinearIdx,
                                 -scaleTwoDt * grid.viscosity(i-1,j));

                    output.addTo(currLinearIdxU,
                                 currLinearIdxU,
                                 scaleTwoDt * grid.viscosity(i-1,j));

                    diag += scaleTwoDt * grid.viscosity(i-1,j);
                    offdiag += std::abs(-scaleTwoDt * grid.viscosity(i-1,j));
                }

                int uIpOneLinearIdx = grid.linearViscosityVelocitySampleIndexU(i+1,j);

                if(uIpOneLinearIdx != -1)
                {
                    output.addTo(currLinearIdxU,
                                 uIpOneLinearIdx,
                                 -scaleTwoDt * grid.viscosity(i,j));

                    output.addTo(currLinearIdxU,
                                 currLinearIdxU,
                                 scaleTwoDt * grid.viscosity(i,j));
                    diag += scaleTwoDt * grid.viscosity(i,j);
                    offdiag += std::abs(-scaleTwoDt * grid.viscosity(i,j));
                }

                int uJmOneLinearIdx = grid.linearViscosityVelocitySampleIndexU(i,j-1);
                int vImOneLinearIdx = grid.linearViscosityVelocitySampleIndexV(i-1,j);

                if(uJmOneLinearIdx != -1
                        && currLinearIdxV != -1
                        && vImOneLinearIdx != -1)
                {
                    float lerpedViscosity = math::lerpCenteredGrid(fi-0.5f,fj-0.5f,grid.viscosityGrid());
                    lerpedViscosity = 100;
                    output.addTo(currLinearIdxU,
                                 uJmOneLinearIdx,
                                 -scaleTwoDx * lerpedViscosity);

                    output.addTo(currLinearIdxU,
                                 vBaseIndex + currLinearIdxV,
                                 scaleTwoDx * lerpedViscosity);

                    output.addTo(currLinearIdxU,
                                 vBaseIndex + vImOneLinearIdx,
                                 -scaleTwoDx * lerpedViscosity);

                    output.addTo(currLinearIdxU,currLinearIdxU,scaleTwoDx * lerpedViscosity);

                    diag += scaleTwoDx * lerpedViscosity;
                    offdiag += std::abs(-scaleTwoDx * lerpedViscosity) * 3;
                }

                int uJpOneLinearIdx = grid.linearViscosityVelocitySampleIndexU(i,j+1);
                int vJpOneLinearIdx = grid.linearViscosityVelocitySampleIndexV(i,j+1);
                int vImOneJpOneLinearIdx = grid.linearViscosityVelocitySampleIndexV(i-1,j+1);

                if(uJpOneLinearIdx != -1
                        && vJpOneLinearIdx != -1
                        && vImOneJpOneLinearIdx != -1)
                {
                    float lerpedViscosity = math::lerpCenteredGrid(fi-0.5f,fj+0.5f,grid.viscosityGrid());
                    lerpedViscosity = 100;
                    output.addTo(currLinearIdxU,
                                 uJpOneLinearIdx,
                                 -scaleTwoDx * lerpedViscosity);

                    output.addTo(currLinearIdxU,
                                 vBaseIndex + vJpOneLinearIdx,
                                 -scaleTwoDx * lerpedViscosity);

                    output.addTo(currLinearIdxU,
                                 vBaseIndex + vImOneJpOneLinearIdx,
                                 scaleTwoDx * lerpedViscosity);

                    output.addTo(currLinearIdxU,currLinearIdxU,scaleTwoDx * lerpedViscosity);

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

            //V component
            if(currLinearIdxV != -1)
            {
                float fi = static_cast<float>(i);
                float fj = static_cast<float>(j);
                int vBaseIndex = validUSamples;
                output.addTo(vBaseIndex + currLinearIdxV,vBaseIndex + currLinearIdxV,SimSettings::density());
                diag += SimSettings::density();

                int vJmOneLinearIdx = grid.linearViscosityVelocitySampleIndexV(i,j-1);

                if(vJmOneLinearIdx != -1)
                {
                    output.addTo(vBaseIndex + currLinearIdxV,
                                 vBaseIndex + vJmOneLinearIdx,
                                 -scaleTwoDt * grid.viscosity(i,j-1));

                    output.addTo(vBaseIndex + currLinearIdxV,
                                 vBaseIndex + currLinearIdxV,
                                 scaleTwoDt * grid.viscosity(i,j-1));
                    diag += scaleTwoDt * grid.viscosity(i,j-1);
                    offdiag += std::abs(scaleTwoDt * grid.viscosity(i,j-1));
                }

                int vJpOneLinearIdx = grid.linearViscosityVelocitySampleIndexV(i,j+1);

                if(grid.vVelocityInside(i,j+1))
                {
                    output.addTo(vBaseIndex + currLinearIdxV,
                                 vBaseIndex + vJpOneLinearIdx,
                                 -scaleTwoDt * grid.viscosity(i,j));

                    output.addTo(vBaseIndex + currLinearIdxV,
                                 vBaseIndex + currLinearIdxV,
                                 scaleTwoDt * grid.viscosity(i,j));

                    diag += scaleTwoDt * grid.viscosity(i,j);
                    offdiag += std::abs(scaleTwoDt * grid.viscosity(i,j));
                }

                int uJmOneLinearIdx = grid.linearViscosityVelocitySampleIndexU(i,j-1);
                int vImOneLinearIdx = grid.linearViscosityVelocitySampleIndexV(i-1,j);

                if(currLinearIdxU != -1
                        && uJmOneLinearIdx != -1
                        && vImOneLinearIdx != -1)
                {
                    float lerpedViscosity = math::lerpCenteredGrid(fi-0.5f,fj-0.5f,grid.viscosityGrid());
                    lerpedViscosity = 100;

                    output.addTo(vBaseIndex + currLinearIdxV,
                                 currLinearIdxU,
                                 scaleTwoDx * lerpedViscosity);

                    output.addTo(vBaseIndex + currLinearIdxV,
                                 uJmOneLinearIdx,
                                 -scaleTwoDx * lerpedViscosity);

                    output.addTo(vBaseIndex + currLinearIdxV,
                                 vBaseIndex + vImOneLinearIdx,
                                 -scaleTwoDx * lerpedViscosity);

                    output.addTo(vBaseIndex + currLinearIdxV,
                                 vBaseIndex + currLinearIdxV,
                                 scaleTwoDx * lerpedViscosity);

                    diag += scaleTwoDx * lerpedViscosity;
                    offdiag += std::abs(scaleTwoDx * lerpedViscosity)*3;
                }

                int uIpOneLinearIdx = grid.linearViscosityVelocitySampleIndexU(i+1,j);
                int uIpOneJmOneLinearIdx = grid.linearViscosityVelocitySampleIndexU(i+1,j-1);
                int vIpOneLinearIdx = grid.linearViscosityVelocitySampleIndexV(i+1,j);

                if(uIpOneLinearIdx != -1
                        && uIpOneJmOneLinearIdx != -1
                        && vIpOneLinearIdx != -1)
                {
                    float lerpedViscosity = math::lerpCenteredGrid(fi+0.5f,fj-0.5f,grid.viscosityGrid());
                    lerpedViscosity = 100;

                    output.addTo(vBaseIndex + currLinearIdxV,
                                 uIpOneLinearIdx,
                                 -scaleTwoDx * lerpedViscosity);

                    output.addTo(vBaseIndex + currLinearIdxV,
                                 uIpOneJmOneLinearIdx,
                                 scaleTwoDx * lerpedViscosity);

                    output.addTo(vBaseIndex + currLinearIdxV,
                                 vBaseIndex + vIpOneLinearIdx,
                                 -scaleTwoDx * lerpedViscosity);

                    output.addTo(vBaseIndex + currLinearIdxV,
                                 vBaseIndex + currLinearIdxV,
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
