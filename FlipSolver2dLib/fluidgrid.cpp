#include "fluidgrid.h"

#include <cmath>

#include "fluidcell.h"
#include "mathfuncs.h"
#include "simsettings.h"

MACFluidGrid::MACFluidGrid(int sizeI, int sizeJ) :
    LinearIndexable2d(sizeI, sizeJ),
    m_materialGrid(sizeI,sizeJ,FluidCellMaterial::EMPTY),
    m_velocityGridU(sizeI + 1,sizeJ, 0.f),
    m_velocityGridV(sizeI, sizeJ + 1, 0.f),
    m_sdf(sizeI,sizeJ,0.f),
    m_knownFlagsU(sizeI + 1, sizeJ, false),
    m_knownFlagsV(sizeI, sizeJ + 1, false),
    m_knownCenteredParams(sizeI,sizeJ, false),
    m_viscosityGrid(sizeI,sizeJ,0.f),
    m_emitterId(sizeI, sizeJ, -1),
    m_solidId(sizeI, sizeJ,-1),
    m_temperature(sizeI, sizeJ, SimSettings::ambientTemp()),
    m_smokeConcentration(sizeI, sizeJ, 0.f),
    m_validUVelocitySampleCount(0),
    m_validVVelocitySampleCount(0),
    m_fluidCellCount(0)
{
}

void MACFluidGrid::fill(FluidCellMaterial m, float velocityU, float velocityV, bool knownFlagU, bool knownFlagV)
{
    if(fluidTest(m))
    {
        m_fluidCellCount = m_materialGrid.sizeI() * m_materialGrid.sizeJ();
    }
    else
    {
        m_fluidCellCount = 0;
    }
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
            if(fluidTest(oldMat) && oldMat != value)
            {
                m_fluidCellCount--;
            }
            else if(!fluidTest(oldMat) && fluidTest(value))
            {
                m_fluidCellCount++;
            }
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

void MACFluidGrid::fillVelocityURect(float value, int topLeftX, int topLeftY, int bottomRightX, int bottomRightY)
{
    m_velocityGridU.fillRect(value,topLeftX,topLeftY,bottomRightX, bottomRightY);
}

void MACFluidGrid::fillVelocityURect(float value, Index2d topLeft, Index2d bottomRight)
{
    m_velocityGridU.fillRect(value,topLeft,bottomRight);
}

void MACFluidGrid::fillVelocityV(float value)
{
    m_velocityGridV.fill(value);
}

void MACFluidGrid::fillVelocityVRect(float value, int topLeftX, int topLeftY, int bottomRightX, int bottomRightY)
{
    m_velocityGridV.fillRect(value,topLeftX,topLeftY,bottomRightX, bottomRightY);
}

void MACFluidGrid::fillVelocityVRect(float value, Index2d topLeft, Index2d bottomRight)
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
    if(!inBounds(index)) return FluidCellMaterial::EMPTY;
    return m_materialGrid.at(index);
}

FluidCellMaterial MACFluidGrid::getMaterial(int i, int j) const
{
    if(!inBounds(i,j)) return FluidCellMaterial::EMPTY;
    return m_materialGrid.at(i,j);
}

void MACFluidGrid::setMaterial(Index2d index, FluidCellMaterial m)
{
    FluidCellMaterial oldMat = getMaterial(index);
    if(fluidTest(oldMat) && !fluidTest(m))
        m_fluidCellCount--;
    else if(!fluidTest(oldMat) && fluidTest(m))
        m_fluidCellCount++;
    m_materialGrid.at(index) = m;
}

void MACFluidGrid::setMaterial(int i, int j, FluidCellMaterial m)
{
    setMaterial(Index2d(i,j),m);
}

bool MACFluidGrid::isFluid(int i, int j)
{
    if(!inBounds(i,j)) return false;
    return fluidTest(m_materialGrid.at(i,j));
}

bool MACFluidGrid::isStrictFluid(int i, int j)
{
    if(!inBounds(i,j)) return false;
    return strictFluidTest(m_materialGrid.at(i,j));
}

bool MACFluidGrid::isSolid(int i, int j)
{
    if(!inBounds(i,j)) return false;
    return solidTest(m_materialGrid.at(i,j));
}

bool MACFluidGrid::isEmpty(int i, int j)
{
    if(!inBounds(i,j)) return true;
    return emptyTest(m_materialGrid.at(i,j));
}

bool MACFluidGrid::isSource(int i, int j)
{
    if(!inBounds(i,j)) return false;
    return sourceTest(m_materialGrid.at(i,j));
}

bool MACFluidGrid::isSink(int i, int j)
{
    if(!inBounds(i,j)) return true;
    return sinkTest(m_materialGrid.at(i,j));
}

bool MACFluidGrid::uVelocitySampleInside(int i, int j)
{
    return (!isEmpty(i,j) && !isEmpty(i-1, j));
}

bool MACFluidGrid::vVelocitySampleInside(int i, int j)
{
    return (!isEmpty(i,j) && !isEmpty(i, j-1));
}

bool MACFluidGrid::uSampleAffectedBySolid(int i, int j)
{
    return (isSolid(i,j+1) || isSolid(i-1, j+1) ||
            isSolid(i,j) || isSolid(i-1, j) ||
            isSolid(i,j-1) || isSolid(i-1, j-1));
}

bool MACFluidGrid::vSampleAffectedBySolid(int i, int j)
{
    return (isSolid(i+1,j) || isSolid(i+1, j-1) ||
            isSolid(i,j) || isSolid(i, j-1) ||
            isSolid(i-1,j) || isSolid(i-1, j-1));
}

VelocitySampleState MACFluidGrid::uVelocitySampleState(int i, int j)
{
    char currentMat = m_materialGrid.at(i,j) >> 4;
    char neighborMat = m_materialGrid.at(i-1,j) >> 4;
    char result = 0;

    result |= currentMat;
    result |= neighborMat;
    if ((result & (result - 1)) != 0)
    {
        result &= ~(0b1000);
    }

    return static_cast<VelocitySampleState>(result);
}

VelocitySampleState MACFluidGrid::vVelocitySampleState(int i, int j)
{
    char currentMat = m_materialGrid.at(i,j) >> 4;
    char neighborMat = m_materialGrid.at(i-1,j) >> 4;
    char result = 0;

    result |= currentMat;
    result |= neighborMat;
    if ((result & (result - 1)) == 0)
    {
        result &= 0b1000;
    }

    return static_cast<VelocitySampleState>(result);
}

void MACFluidGrid::updateValidULinearMapping()
{
    m_uVelocitySamplesMap.clear();
    int linearIdx = 0;
    for(int i = 0; i < m_sizeI; i++)
    {
        for(int j = 0; j < m_sizeJ; j++)
        {
            if(uVelocitySampleInside(i,j))
            {
                m_uVelocitySamplesMap.insert(std::make_pair(std::pair(i,j),linearIdx));
                linearIdx++;
            }
        }
    }

    m_validUVelocitySampleCount = linearIdx;
}

void MACFluidGrid::updateValidVLinearMapping()
{
    m_vVelocitySamplesMap.clear();
    int linearIdx = 0;
    for(int i = 0; i < m_sizeI; i++)
    {
        for(int j = 0; j < m_sizeJ; j++)
        {
            if(vVelocitySampleInside(i,j))
            {
                m_vVelocitySamplesMap.insert(std::make_pair(std::pair(i,j),linearIdx));
                linearIdx++;
            }
        }
    }

    m_validVVelocitySampleCount = linearIdx;
}

void MACFluidGrid::setU(Index2d index, float value, bool knownStatus)
{
    m_velocityGridU.at(index) = value;
    m_knownFlagsU.at(index) = knownStatus;
}

void MACFluidGrid::setU(Index2d index, float value)
{
    m_velocityGridU.at(index) = value;
}

void MACFluidGrid::setU(int i, int j, float value, bool knownStatus)
{
    m_velocityGridU.at(i,j) = value;
    m_knownFlagsU.at(i,j) = knownStatus;
}

void MACFluidGrid::setU(int i, int j, float value)
{
    m_velocityGridU.at(i,j) = value;
}

void MACFluidGrid::setV(Index2d index, float value, bool knownStatus)
{
    m_velocityGridV.at(index) = value;
    m_knownFlagsV.at(index) = knownStatus;
}

void MACFluidGrid::setV(Index2d index, float value)
{
    m_velocityGridV.at(index) = value;
}

void MACFluidGrid::setV(int i, int j, float value, bool knownStatus)
{
    m_velocityGridV.at(i,j) = value;
    m_knownFlagsV.at(i,j) = knownStatus;
}

void MACFluidGrid::setV(int i, int j, float value)
{
    m_velocityGridV.at(i,j) = value;
}

float MACFluidGrid::getU(Index2d index) const
{
    return m_velocityGridU.at(index);
}

float MACFluidGrid::getU(int i, int j) const
{
    return m_velocityGridU.at(i,j);
}

float MACFluidGrid::getV(Index2d index) const
{
    return m_velocityGridV.at(index);
}

float MACFluidGrid::getV(int i, int j) const
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

Vertex MACFluidGrid::velocityAt(float i, float j)
{
    return Vertex(math::lerpUGrid(i, j, m_velocityGridU),math::lerpVGrid(i, j, m_velocityGridV));
}

Vertex MACFluidGrid::velocityAt(Vertex position)
{
    return velocityAt(position.x(),position.y());
}

float MACFluidGrid::viscosityAt(Vertex position)
{
    return math::lerpCenteredGrid(position.x(), position.y(), m_viscosityGrid);
}

Grid2d<FluidCellMaterial> &MACFluidGrid::materialGrid()
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

Grid2d<bool> &MACFluidGrid::knownFlagsCenteredParams()
{
    return m_knownCenteredParams;
}

Grid2d<float> &MACFluidGrid::viscosityGrid()
{
    return m_viscosityGrid;
}

Grid2d<float> &MACFluidGrid::sdfGrid()
{
    return m_sdf;
}

Grid2d<float> &MACFluidGrid::temperatureGrid()
{
    return m_temperature;
}

Grid2d<float> &MACFluidGrid::smokeConcentrationGrid()
{
    return m_smokeConcentration;
}

float MACFluidGrid::sdf(int i, int j)
{
    return m_sdf.at(i,j);
}

float MACFluidGrid::sdfAt(float i, float j)
{
    return math::lerpCenteredGrid(i,j,m_sdf);
}

float MACFluidGrid::sdf(Index2d index)
{
    return m_sdf.at(index);
}

void MACFluidGrid::setSdf(int i, int j, float value)
{
    m_sdf.at(i,j) = value;
}

void MACFluidGrid::setSdf(Index2d index, float value)
{
    m_sdf.at(index) = value;
}

float MACFluidGrid::viscosity(int i, int j)
{
    return m_viscosityGrid.at(i,j);
}

float MACFluidGrid::viscosityAt(float i, float j)
{
    return math::lerpCenteredGrid(i,j,m_viscosityGrid);
}

float MACFluidGrid::viscosity(Index2d index)
{
    return m_viscosityGrid.at(index);
}

void MACFluidGrid::setViscosity(int i, int j, float value)
{
    m_viscosityGrid.at(i,j) = value;
}

void MACFluidGrid::setViscosity(Index2d index, float value)
{
    m_viscosityGrid.at(index) = value;
}

float MACFluidGrid::temperature(int i, int j)
{
    return m_temperature.at(i,j);
}

float MACFluidGrid::temperatureAt(float i, float j)
{
    return math::lerpCenteredGrid(i,j,m_temperature);
}

void MACFluidGrid::setTemperature(int i, int j, float value)
{
    m_temperature.setAt(i,j,value);
}

float MACFluidGrid::smokeConcentration(int i, int j)
{
    return m_smokeConcentration.at(i,j);
}

float MACFluidGrid::smokeConcentrationAt(float i, float j)
{
    return math::lerpCenteredGrid(i,j,m_smokeConcentration);
}

void MACFluidGrid::setSmokeConcentration(int i, int j, float value)
{
    m_smokeConcentration.setAt(i,j,value);
}

int MACFluidGrid::closestSolidId(int i, int j)
{
    Vertex surfacePoint = closestSurfacePoint(Vertex(static_cast<float>(i) + 0.5f,
                                              static_cast<float>(j) + 0.5f));
    return m_solidId.at(static_cast<int>(surfacePoint.x()),static_cast<int>(surfacePoint.y()));
}

std::vector<int> MACFluidGrid::nearValidSolidIds(int i, int j)
{
    std::vector<int> output;
    for (int iOffset = -1; iOffset < 2; iOffset++)
    {
        for (int jOffset = -1; jOffset < 2; jOffset++)
        {
            if(m_solidId.inBounds(i + iOffset, j + jOffset))
            {
                int id = m_solidId.at(i + iOffset, j + jOffset);
                if(id != -1)
                {
                    output.push_back(m_solidId.at(i + iOffset, j + jOffset));
                }
            }
        }
    }

    return output;
}

Vertex MACFluidGrid::closestSurfacePoint(Vertex pos)
{
    Vertex closestPoint = pos;
    float sdf = math::lerpCenteredGrid(pos.x(),pos.y(),m_sdf);
    Vertex grad = math::gradCenteredGrid(pos.x(),pos.y(),m_sdf);
    static const int iterationLimit = 100;
    static const int internalIterationLimit = 10;
    for(int i = 0; i < iterationLimit; i++)
    {
        float alpha = 1;
        for(int j = 0; j < internalIterationLimit; j++)
        {
            Vertex q = closestPoint - alpha*sdf*grad;
            if(std::abs(math::lerpCenteredGrid(q.x(),q.y(),m_sdf)) < std::abs(sdf))
            {
                closestPoint = q;
                sdf = math::lerpCenteredGrid(q.x(),q.y(),m_sdf);
                grad = math::gradCenteredGrid(q.x(),q.y(),m_sdf);
                if(std::abs(sdf) < 1e-5f)
                {
                    return closestPoint;
                }
            }
            else
            {
                alpha *= 0.7f;
            }
        }
    }
    return closestPoint;
}

void MACFluidGrid::updateLinearFluidViscosityMapping()
{
    updateValidULinearMapping();
    updateValidVLinearMapping();
}

int MACFluidGrid::linearViscosityVelocitySampleIndexU(int i, int j)
{
    std::unordered_map<std::pair<int,int>,int>::iterator iter = m_uVelocitySamplesMap.find(std::pair<int,int>(i,j));
    if(iter != m_uVelocitySamplesMap.end())
    {
        return iter->second;
    }
    else
    {
        return -1;
    }
}

int MACFluidGrid::linearViscosityVelocitySampleIndexV(int i, int j)
{
    std::unordered_map<std::pair<int,int>,int>::iterator iter = m_vVelocitySamplesMap.find(std::pair<int,int>(i,j));
    if(iter != m_vVelocitySamplesMap.end())
    {
        return iter->second;
    }
    else
    {
        return -1;
    }
}

int MACFluidGrid::validVelocitySampleCountU()
{
    return m_validUVelocitySampleCount;
}

int MACFluidGrid::validVelocitySampleCountV()
{
    return m_validVVelocitySampleCount;
}

void MACFluidGrid::setEmitterId(int i, int j, int id)
{
    m_emitterId.at(i,j) = id;
}

int MACFluidGrid::emitterId(int i, int j)
{
    return m_emitterId.at(i,j);
}

void MACFluidGrid::setSolidId(int i, int j, int solidId)
{
    m_solidId.at(i,j) = solidId;
}

int MACFluidGrid::solidId(int i, int j)
{
    return m_solidId.at(i,j);
}
