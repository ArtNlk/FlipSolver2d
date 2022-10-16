#include "fluidgrid.h"

#include <cmath>

#include "fluidcell.h"
#include "mathfuncs.h"
#include "simsettings.h"

MACFluidGrid::MACFluidGrid(int sizeI, int sizeJ) :
    LinearIndexable2d(sizeI, sizeJ),
    m_materialGrid(sizeI,sizeJ,FluidMaterial::EMPTY),
    m_fluidVelocityGrid(sizeI, sizeJ),
    m_airVelocityGrid(sizeI, sizeJ),
    m_sdf(sizeI,sizeJ,0.f),
    m_knownCenteredParams(sizeI,sizeJ, false),
    m_viscosityGrid(sizeI,sizeJ,0.f),
    m_emitterId(sizeI, sizeJ, -1),
    m_solidId(sizeI, sizeJ,-1),
    m_temperature(sizeI, sizeJ, SimSettings::ambientTemp()),
    m_smokeConcentration(sizeI, sizeJ, 0.f),
    m_divergenceControl(sizeI,sizeJ, 0.f),
    m_particleCounts(sizeI, sizeJ),
    m_validUVelocitySampleCount(0),
    m_validVVelocitySampleCount(0),
    m_fluidCellCount(0)
{
}

void MACFluidGrid::fill(FluidMaterial m, float velocityU, float velocityV, bool knownFlagU, bool knownFlagV)
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

void MACFluidGrid::fillMaterial(FluidMaterial m)
{
    m_materialGrid.fill(m);
}

void MACFluidGrid::fillMaterialRect(FluidMaterial value, int topLeftX, int topLeftY, int bottomRightX, int bottomRightY)
{
    topLeftX = std::max(0,topLeftX);
    topLeftY = std::max(0,topLeftY);
    bottomRightX = std::min(m_sizeI - 1,bottomRightX);
    bottomRightY = std::min(m_sizeJ - 1,bottomRightY);
    for(int i = topLeftX; i <= bottomRightX; i++)
    {
        for(int j = topLeftY; j <= bottomRightY; j++)
        {
            FluidMaterial oldMat = m_materialGrid.at(i,j);
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

void MACFluidGrid::fillMaterialRect(FluidMaterial value, Index2d topLeft, Index2d bottomRight)
{
    fillMaterialRect(value,topLeft.m_i,topLeft.m_j,bottomRight.m_i,bottomRight.m_j);
}

void MACFluidGrid::fillVelocityU(float value)
{
    m_fluidVelocityGrid.velocityGridU().fill(value);
}

void MACFluidGrid::fillVelocityURect(float value, int topLeftX, int topLeftY, int bottomRightX, int bottomRightY)
{
    m_fluidVelocityGrid.velocityGridU().fillRect(value,topLeftX,topLeftY,bottomRightX, bottomRightY);
}

void MACFluidGrid::fillVelocityURect(float value, Index2d topLeft, Index2d bottomRight)
{
    m_fluidVelocityGrid.velocityGridU().fillRect(value,topLeft,bottomRight);
}

void MACFluidGrid::fillVelocityV(float value)
{
    m_fluidVelocityGrid.velocityGridV().fill(value);
}

void MACFluidGrid::fillVelocityVRect(float value, int topLeftX, int topLeftY, int bottomRightX, int bottomRightY)
{
    m_fluidVelocityGrid.velocityGridV().fillRect(value,topLeftX,topLeftY,bottomRightX, bottomRightY);
}

void MACFluidGrid::fillVelocityVRect(float value, Index2d topLeft, Index2d bottomRight)
{
    m_fluidVelocityGrid.velocityGridV().fillRect(value,topLeft,bottomRight);
}

void MACFluidGrid::fillKnownFlagsU(bool value)
{
    m_fluidVelocityGrid.uSampleValidityGrid().fill(value);
}

void MACFluidGrid::fillKnownFlagURect(bool value, int topLeftX, int topLeftY, int bottomRightX, int bottomRightY)
{
    m_fluidVelocityGrid.uSampleValidityGrid().fillRect(value,topLeftX,topLeftY,bottomRightX, bottomRightY);
}

void MACFluidGrid::fillKnownFlagURect(bool value, Index2d topLeft, Index2d bottomRight)
{
    m_fluidVelocityGrid.uSampleValidityGrid().fillRect(value,topLeft,bottomRight);
}

void MACFluidGrid::fillKnownFlagsV(bool value)
{
    m_fluidVelocityGrid.vSampleValidityGrid().fill(value);
}

void MACFluidGrid::fillKnownFlagVRect(bool value, int topLeftX, int topLeftY, int bottomRightX, int bottomRightY)
{
    m_fluidVelocityGrid.vSampleValidityGrid().fillRect(value,topLeftX,topLeftY,bottomRightX, bottomRightY);
}

void MACFluidGrid::fillKnownFlagVRect(bool value, Index2d topLeft, Index2d bottomRight)
{
    m_fluidVelocityGrid.vSampleValidityGrid().fillRect(value,topLeft,bottomRight);
}

FluidMaterial MACFluidGrid::getMaterial(Index2d index) const
{
    if(!inBounds(index)) return FluidMaterial::EMPTY;
    return m_materialGrid.at(index);
}

FluidMaterial MACFluidGrid::getMaterial(int i, int j) const
{
    if(!inBounds(i,j)) return FluidMaterial::EMPTY;
    return m_materialGrid.at(i,j);
}

void MACFluidGrid::setMaterial(Index2d index, FluidMaterial m)
{
    FluidMaterial oldMat = getMaterial(index);
    if(fluidTest(oldMat) && !fluidTest(m))
        m_fluidCellCount--;
    else if(!fluidTest(oldMat) && fluidTest(m))
        m_fluidCellCount++;
    m_materialGrid.at(index) = m;
}

void MACFluidGrid::setMaterial(int i, int j, FluidMaterial m)
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

void MACFluidGrid::setFluidU(Index2d index, float value, bool knownStatus)
{
    m_fluidVelocityGrid.velocityGridU().at(index) = value;
    m_fluidVelocityGrid.uSampleValidityGrid().at(index) = knownStatus;
}

void MACFluidGrid::setFluidU(Index2d index, float value)
{
    m_fluidVelocityGrid.velocityGridU().at(index) = value;
}

void MACFluidGrid::setFluidU(int i, int j, float value, bool knownStatus)
{
    m_fluidVelocityGrid.velocityGridU().at(i,j) = value;
    m_fluidVelocityGrid.uSampleValidityGrid().at(i,j) = knownStatus;
}

void MACFluidGrid::setFluidU(int i, int j, float value)
{
    m_fluidVelocityGrid.velocityGridU().at(i,j) = value;
}

void MACFluidGrid::setFluidV(Index2d index, float value, bool knownStatus)
{
    m_fluidVelocityGrid.velocityGridV().at(index) = value;
    m_fluidVelocityGrid.vSampleValidityGrid().at(index) = knownStatus;
}

void MACFluidGrid::setFluidV(Index2d index, float value)
{
    m_fluidVelocityGrid.velocityGridV().at(index) = value;
}

void MACFluidGrid::setFluidV(int i, int j, float value, bool knownStatus)
{
    m_fluidVelocityGrid.velocityGridV().at(i,j) = value;
    m_fluidVelocityGrid.vSampleValidityGrid().at(i,j) = knownStatus;
}

void MACFluidGrid::setFluidV(int i, int j, float value)
{
    m_fluidVelocityGrid.velocityGridV().at(i,j) = value;
}

void MACFluidGrid::setAirU(int i, int j, float value)
{
    m_airVelocityGrid.velocityGridU().at(i,j) = value;
}

void MACFluidGrid::setAirV(int i, int j, float value)
{
    m_airVelocityGrid.velocityGridV().at(i,j) = value;
}

void MACFluidGrid::getSize(int &sizeI, int &sizeJ) const
{
    sizeI = m_sizeI;
    sizeJ = m_sizeJ;
}

float MACFluidGrid::getFluidU(int i, int j)
{
    return m_fluidVelocityGrid.velocityGridU().at(i,j);
}

float MACFluidGrid::getFluidV(int i, int j)
{
    return m_fluidVelocityGrid.velocityGridV().at(i,j);
}

float MACFluidGrid::getAirU(int i, int j)
{
    return m_airVelocityGrid.velocityGridU().at(i,j);
}

float MACFluidGrid::getAirV(int i, int j)
{
    return m_airVelocityGrid.velocityGridV().at(i,j);
}

int MACFluidGrid::cellCount() const
{
    return m_sizeI * m_sizeJ;
}

int MACFluidGrid::fluidCellCount() const
{
    return m_fluidCellCount;
}

Vertex MACFluidGrid::fluidVelocityAt(float i, float j)
{
    return m_fluidVelocityGrid.velocityAt(i,j);
}

Vertex MACFluidGrid::fluidVelocityAt(Vertex position)
{
    return m_fluidVelocityGrid.velocityAt(position.x(),position.y());
}

float MACFluidGrid::viscosityAt(Vertex position)
{
    return simmath::lerpCenteredGrid(position.x(), position.y(), m_viscosityGrid);
}

Grid2d<FluidMaterial> &MACFluidGrid::materialGrid()
{
    return m_materialGrid;
}

StaggeredVelocityGrid &MACFluidGrid::fluidVelocityGrid()
{
    return m_fluidVelocityGrid;
}

StaggeredVelocityGrid &MACFluidGrid::airVelocityGrid()
{
    return m_airVelocityGrid;
}

Grid2d<float> &MACFluidGrid::velocityGridU()
{
    return m_fluidVelocityGrid.velocityGridU();
}

Grid2d<float> &MACFluidGrid::velocityGridV()
{
    return m_fluidVelocityGrid.velocityGridV();
}

Grid2d<bool> &MACFluidGrid::knownFlagsGridU()
{
    return m_fluidVelocityGrid.uSampleValidityGrid();
}

Grid2d<bool> &MACFluidGrid::knownFlagsGridV()
{
    return m_fluidVelocityGrid.vSampleValidityGrid();
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

Grid2d<float> &MACFluidGrid::divergenceControlGrid()
{
    return m_divergenceControl;
}

Grid2d<int> &MACFluidGrid::particleCountGrid()
{
    return m_particleCounts;
}

float MACFluidGrid::sdf(int i, int j)
{
    return m_sdf.at(i,j);
}

float MACFluidGrid::sdfAt(float i, float j)
{
    return simmath::lerpCenteredGrid(i,j,m_sdf);
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
    return simmath::lerpCenteredGrid(i,j,m_viscosityGrid);
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
    return simmath::lerpCenteredGrid(i,j,m_temperature);
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
    return simmath::lerpCenteredGrid(i,j,m_smokeConcentration);
}

void MACFluidGrid::setSmokeConcentration(int i, int j, float value)
{
    m_smokeConcentration.setAt(i,j,value);
}

float& MACFluidGrid::divergenceControl(int i, int j)
{
    return m_divergenceControl.at(i,j);
}

void MACFluidGrid::setDivergeceControl(int i, int j, float value)
{
    m_divergenceControl.at(i,j) = value;
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
    float sdf = simmath::lerpCenteredGrid(pos.x(),pos.y(),m_sdf);
    Vertex grad = simmath::gradCenteredGrid(pos.x(),pos.y(),m_sdf);
    static const int iterationLimit = 100;
    static const int internalIterationLimit = 10;
    for(int i = 0; i < iterationLimit; i++)
    {
        float alpha = 1;
        for(int j = 0; j < internalIterationLimit; j++)
        {
            Vertex q = closestPoint - alpha*sdf*grad;
            if(std::abs(simmath::lerpCenteredGrid(q.x(),q.y(),m_sdf)) < std::abs(sdf))
            {
                closestPoint = q;
                sdf = simmath::lerpCenteredGrid(q.x(),q.y(),m_sdf);
                grad = simmath::gradCenteredGrid(q.x(),q.y(),m_sdf);
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
