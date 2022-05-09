#include "fluidgrid.h"

MACFluidGrid::MACFluidGrid(int sizeI, int sizeJ) :
    LinearIndexable2d(sizeI, sizeJ),
    m_materialGrid(sizeI,sizeJ,FluidCellMaterial::EMPTY),
    m_velocityGridU(sizeI + 1,sizeJ, 0.f),
    m_velocityGridV(sizeI, sizeJ + 1, 0.f),
    m_knownFlagsU(sizeI + 1, sizeJ, false),
    m_knownFlagsV(sizeI, sizeJ + 1, false),
    m_fluidCellCount(0)
{
}

void MACFluidGrid::fill(FluidCellMaterial m, float velocityU, float velocityV, bool knownFlagU, bool knownFlagV)
{
    if(m == FluidCellMaterial::FLUID)
        m_fluidCellCount = m_materialGrid.sizeI() * m_materialGrid.sizeJ();
    else
        m_fluidCellCount = 0;
    fillMaterial(m);
    fillVelocityU(velocityU);
    fillVelocityV(velocityV);
    fillKnownFlagsU(knownFlagU);
    fillKnownFlagsV(knownFlagV);
}

void MACFluidGrid::fillMaterial(FluidCellMaterial m)
{
    m_materialGrid.fill(m);
}

void MACFluidGrid::fillMaterialRect(FluidCellMaterial value, int topLeftX, int topLeftY, int bottomRightX, int bottomRightY)
{
    topLeftX = std::max(0,topLeftX);
    topLeftY = std::max(0,topLeftY);
    bottomRightX = std::min(m_sizeI - 1,bottomRightX);
    bottomRightY = std::min(m_sizeJ - 1,bottomRightY);
    for(int i = topLeftX; i <= bottomRightX; i++)
    {
        for(int j = topLeftY; j <= bottomRightY; j++)
        {
            FluidCellMaterial oldMat = m_materialGrid.at(i,j);
            if(oldMat == FluidCellMaterial::FLUID && oldMat != value)
                m_fluidCellCount--;
            else if(oldMat != FluidCellMaterial::FLUID && value == FluidCellMaterial::FLUID)
                m_fluidCellCount++;
            m_materialGrid.setAt(i,j,value);
        }
    }
}

void MACFluidGrid::fillMaterialRect(FluidCellMaterial value, Index2d topLeft, Index2d bottomRight)
{
    fillMaterialRect(value,topLeft.m_i,topLeft.m_j,bottomRight.m_i,bottomRight.m_j);
}

void MACFluidGrid::fillVelocityU(float value)
{
    m_velocityGridU.fill(value);
}

void MACFluidGrid::fillVelocityURect(double value, int topLeftX, int topLeftY, int bottomRightX, int bottomRightY)
{
    m_velocityGridU.fillRect(value,topLeftX,topLeftY,bottomRightX, bottomRightY);
}

void MACFluidGrid::fillVelocityURect(double value, Index2d topLeft, Index2d bottomRight)
{
    m_velocityGridU.fillRect(value,topLeft,bottomRight);
}

void MACFluidGrid::fillVelocityV(float value)
{
    m_velocityGridV.fill(value);
}

void MACFluidGrid::fillVelocityVRect(double value, int topLeftX, int topLeftY, int bottomRightX, int bottomRightY)
{
    m_velocityGridV.fillRect(value,topLeftX,topLeftY,bottomRightX, bottomRightY);
}

void MACFluidGrid::fillVelocityVRect(double value, Index2d topLeft, Index2d bottomRight)
{
    m_velocityGridV.fillRect(value,topLeft,bottomRight);
}

void MACFluidGrid::fillKnownFlagsU(bool value)
{
    m_knownFlagsU.fill(value);
}

void MACFluidGrid::fillKnownFlagURect(bool value, int topLeftX, int topLeftY, int bottomRightX, int bottomRightY)
{
    m_knownFlagsU.fillRect(value,topLeftX,topLeftY,bottomRightX, bottomRightY);
}

void MACFluidGrid::fillKnownFlagURect(bool value, Index2d topLeft, Index2d bottomRight)
{
    m_knownFlagsU.fillRect(value,topLeft,bottomRight);
}

void MACFluidGrid::fillKnownFlagsV(bool value)
{
    m_knownFlagsV.fill(value);
}

void MACFluidGrid::fillKnownFlagVRect(bool value, int topLeftX, int topLeftY, int bottomRightX, int bottomRightY)
{
    m_knownFlagsV.fillRect(value,topLeftX,topLeftY,bottomRightX, bottomRightY);
}

void MACFluidGrid::fillKnownFlagVRect(bool value, Index2d topLeft, Index2d bottomRight)
{
    m_knownFlagsV.fillRect(value,topLeft,bottomRight);
}

FluidCellMaterial MACFluidGrid::getMaterial(Index2d index) const
{
    if(!inBounds(index)) return FluidCellMaterial::SOLID;
    return m_materialGrid.at(index);
}

FluidCellMaterial MACFluidGrid::getMaterial(int i, int j) const
{
    if(!inBounds(i,j)) return FluidCellMaterial::SOLID;
    return m_materialGrid.at(i,j);
}

void MACFluidGrid::setMaterial(Index2d index, FluidCellMaterial m)
{
    FluidCellMaterial oldMat = getMaterial(index);
    if(oldMat == FluidCellMaterial::FLUID && oldMat != m)
        m_fluidCellCount--;
    else if(oldMat != FluidCellMaterial::FLUID && m == FluidCellMaterial::FLUID)
        m_fluidCellCount++;
    m_materialGrid.at(index) = m;
}

void MACFluidGrid::setMaterial(int i, int j, FluidCellMaterial m)
{
    FluidCellMaterial oldMat = getMaterial(i,j);
    if(oldMat == FluidCellMaterial::FLUID && oldMat != m)
        m_fluidCellCount--;
    else if(oldMat != FluidCellMaterial::FLUID && m == FluidCellMaterial::FLUID)
        m_fluidCellCount++;
    m_materialGrid.at(i,j) = m;
}

void MACFluidGrid::setU(Index2d index, double value, bool knownStatus)
{
    m_velocityGridU.at(index) = value;
    m_knownFlagsU.at(index) = knownStatus;
}

void MACFluidGrid::setU(Index2d index, double value)
{
    m_velocityGridU.at(index) = value;
}

void MACFluidGrid::setU(int i, int j, double value, bool knownStatus)
{
    m_velocityGridU.at(i,j) = value;
    m_knownFlagsU.at(i,j) = knownStatus;
}

void MACFluidGrid::setU(int i, int j, double value)
{
    m_velocityGridU.at(i,j) = value;
}

void MACFluidGrid::setV(Index2d index, double value, bool knownStatus)
{
    m_velocityGridV.at(index) = value;
    m_knownFlagsV.at(index) = knownStatus;
}

void MACFluidGrid::setV(Index2d index, double value)
{
    m_velocityGridV.at(index) = value;
}

void MACFluidGrid::setV(int i, int j, double value, bool knownStatus)
{
    m_velocityGridV.at(i,j) = value;
    m_knownFlagsV.at(i,j) = knownStatus;
}

void MACFluidGrid::setV(int i, int j, double value)
{
    m_velocityGridV.at(i,j) = value;
}

double MACFluidGrid::getU(Index2d index) const
{
    return m_velocityGridU.at(index);
}

double MACFluidGrid::getU(int i, int j) const
{
    return m_velocityGridU.at(i,j);
}

double MACFluidGrid::getV(Index2d index) const
{
    return m_velocityGridV.at(index);
}

double MACFluidGrid::getV(int i, int j) const
{
    return m_velocityGridV.at(i,j);
}

void MACFluidGrid::getSize(int &sizeI, int &sizeJ) const
{
    sizeI = m_sizeI;
    sizeJ = m_sizeJ;
}

int MACFluidGrid::cellCount() const
{
    return m_sizeI * m_sizeJ;
}

int MACFluidGrid::fluidCellCount() const
{
    return m_fluidCellCount;
}

const Grid2d<FluidCellMaterial> &MACFluidGrid::materialGrid()
{
    return m_materialGrid;
}

Grid2d<float> &MACFluidGrid::velocityGridU()
{
    return m_velocityGridU;
}

Grid2d<float> &MACFluidGrid::velocityGridV()
{
    return m_velocityGridV;
}

Grid2d<bool> &MACFluidGrid::knownFlagsGridU()
{
    return m_knownFlagsU;
}

Grid2d<bool> &MACFluidGrid::knownFlagsGridV()
{
    return m_knownFlagsV;
}

void MACFluidGrid::updateLinearToFluidMapping()
{
    m_linearToFluidCellIndexMap.clear();
    int fluidCellIndex = 0;
    for(int i = 0; i < m_sizeI; i++)
    {
        for(int j = 0; j < m_sizeJ; j++)
        {
            if(getMaterial(i,j) == FluidCellMaterial::FLUID)
            {
                m_linearToFluidCellIndexMap[linearIndex(i,j)] = fluidCellIndex;
                fluidCellIndex++;
            }
        }
    }
}

int MACFluidGrid::linearFluidIndex(int i, int j)
{
    std::unordered_map<int,int>::iterator iter = m_linearToFluidCellIndexMap.find(linearIndex(i,j));
    if(iter != m_linearToFluidCellIndexMap.end())
    {
        return iter->second;
    }
    else
    {
        return -1;
    }
}
