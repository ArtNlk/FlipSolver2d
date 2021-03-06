#include "fluidgrid.h"

#include <cmath>

#include "mathfuncs.h"

MACFluidGrid::MACFluidGrid(int sizeI, int sizeJ) :
    LinearIndexable2d(sizeI, sizeJ),
    m_materialGrid(sizeI,sizeJ,FluidCellMaterial::EMPTY),
    m_velocityGridU(sizeI + 1,sizeJ, 0.f),
    m_velocityGridV(sizeI, sizeJ + 1, 0.f),
    m_sdf(sizeI,sizeJ,0.f),
    m_knownFlagsU(sizeI + 1, sizeJ, false),
    m_knownFlagsV(sizeI, sizeJ + 1, false),
    m_viscosity(sizeI,sizeJ,0.f),
    m_fluidCellCount(0)
{
}

void MACFluidGrid::fill(FluidCellMaterial m, float velocityU, float velocityV, bool knownFlagU, bool knownFlagV)
{
    if(fluidTest(m))
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
            if(fluidTest(oldMat) && oldMat != value)
                m_fluidCellCount--;
            else if(!fluidTest(oldMat) && fluidTest(value))
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

void MACFluidGrid::getFlattenedFluidVelocities(std::vector<double> &velocities)
{
    ASSERT(velocities.size() == (m_fluidCellCount * 2));
    for(std::pair<std::pair<int,int>,int> p : m_fluidCellIndexMap)
    {
        velocities[p.second] = m_velocityGridU.at(p.first.first,p.first.second);
        velocities[m_fluidCellCount + p.second] = m_velocityGridV.at(p.first.first,p.first.second);
    }
}

void MACFluidGrid::unflattenFluidVelocities(std::vector<double> &velocities)
{
    ASSERT(velocities.size() == (m_fluidCellCount * 2));
    for(std::pair<std::pair<int,int>,int> p : m_fluidCellIndexMap)
    {
        m_velocityGridU.at(p.first.first,p.first.second) = velocities[p.second];
        m_velocityGridV.at(p.first.first,p.first.second) = velocities[m_fluidCellCount + p.second];
    }
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

Grid2d<float> &MACFluidGrid::viscosityGrid()
{

}

Grid2d<float> &MACFluidGrid::sdfGrid()
{
    return m_sdf;
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
    return m_viscosity.at(i,j);
}

float MACFluidGrid::viscosityAt(float i, float j)
{
    return math::lerpCenteredGrid(i,j,m_viscosity);
}

float MACFluidGrid::viscosity(Index2d index)
{
    return m_viscosity.at(index);
}

void MACFluidGrid::setViscosity(int i, int j, float value)
{
    m_viscosity.at(i,j) = value;
}

void MACFluidGrid::setViscosity(Index2d index, float value)
{
    m_viscosity.at(index) = value;
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
            float temp = std::abs(math::lerpCenteredGrid(q.x(),q.y(),m_sdf));
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

void MACFluidGrid::updateLinearFluidCellMapping()
{
    m_fluidCellIndexMap.clear();
    m_fluidCellCount = 0;
    for(int i = 0; i < m_sizeI; i++)
    {
        for(int j = 0; j < m_sizeJ; j++)
        {
            if(fluidTest(getMaterial(i,j)))
            {
                m_fluidCellIndexMap[std::pair<int,int>(i,j)] = m_fluidCellCount;
                m_fluidCellCount++;
            }
        }
    }
}

int MACFluidGrid::linearFluidIndex(int i, int j)
{
    std::unordered_map<std::pair<int,int>,int>::iterator iter = m_fluidCellIndexMap.find(std::pair<int,int>(i,j));
    if(iter != m_fluidCellIndexMap.end())
    {
        return iter->second;
    }
    else
    {
        return -1;
    }
}
