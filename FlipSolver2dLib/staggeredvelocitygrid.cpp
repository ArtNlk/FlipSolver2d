#include "staggeredvelocitygrid.h"

#include "grid2d.h"
#include "mathfuncs.h"

StaggeredVelocityGrid::StaggeredVelocityGrid(size_t sizeI, size_t sizeJ) :
    LinearIndexable2d(sizeI, sizeJ),
    m_velocityGridU(sizeI + 1, sizeJ, 0.f, OOBStrategy::OOB_EXTEND, 0.f, Vec3(0.5f,0.f)),
    m_velocityGridV(sizeI, sizeJ + 1, 0.f, OOBStrategy::OOB_EXTEND, 0.f, Vec3(0.f, 0.5f)),
    m_uSampleValidity(sizeI + 1, sizeJ, false, OOBStrategy::OOB_CONST, true),
    m_vSampleValidity(sizeI, sizeJ + 1, false, OOBStrategy::OOB_CONST, true)
{

}

Grid2d<float> &StaggeredVelocityGrid::velocityGridU()
{
    return m_velocityGridU;
}

Grid2d<float> &StaggeredVelocityGrid::velocityGridV()
{
    return m_velocityGridV;
}

const Grid2d<float> &StaggeredVelocityGrid::velocityGridU() const
{
    return m_velocityGridU;
}

const Grid2d<float> &StaggeredVelocityGrid::velocityGridV() const
{
    return m_velocityGridV;
}

Grid2d<bool> &StaggeredVelocityGrid::uSampleValidityGrid()
{
    return m_uSampleValidity;
}

Grid2d<bool> &StaggeredVelocityGrid::vSampleValidityGrid()
{
    return m_vSampleValidity;
}

float &StaggeredVelocityGrid::u(ssize_t i, ssize_t j)
{
    return m_velocityGridU.at(i,j);
}

float &StaggeredVelocityGrid::v(ssize_t i, ssize_t j)
{
    return m_velocityGridV.at(i,j);
}

float StaggeredVelocityGrid::getU(ssize_t i, ssize_t j) const
{
    return m_velocityGridU.getAt(i,j);
}

float StaggeredVelocityGrid::getV(ssize_t i, ssize_t j) const
{
    return m_velocityGridV.getAt(i,j);
}

void StaggeredVelocityGrid::setU(ssize_t i, ssize_t j, float u)
{
    m_velocityGridU.setAt(i,j,u);
}

void StaggeredVelocityGrid::setV(ssize_t i, ssize_t j, float v)
{
    m_velocityGridV.setAt(i,j,v);
}

bool StaggeredVelocityGrid::getUValidity(ssize_t i, ssize_t j)
{
    return m_uSampleValidity.getAt(i,j);
}

bool StaggeredVelocityGrid::getVValidity(ssize_t i, ssize_t j)
{
    return m_vSampleValidity.getAt(i,j);
}

void StaggeredVelocityGrid::setUValidity(ssize_t i, ssize_t j, bool uValidity)
{
    m_uSampleValidity.setAt(i,j, uValidity);
}

void StaggeredVelocityGrid::setVValidity(ssize_t i, ssize_t j, bool vValidity)
{
    m_vSampleValidity.setAt(i,j, vValidity);
}

float &StaggeredVelocityGrid::u(Index2d idx)
{
    return u(idx.i, idx.j);
}

float &StaggeredVelocityGrid::v(Index2d idx)
{
    return v(idx.i, idx.j);
}

float StaggeredVelocityGrid::getU(Index2d idx) const
{
    return getU(idx.i, idx.j);
}

float StaggeredVelocityGrid::getV(Index2d idx) const
{
    return getV(idx.i, idx.j);
}

void StaggeredVelocityGrid::setU(Index2d idx, float u)
{
    setU(idx.i, idx.j,u);
}

void StaggeredVelocityGrid::setV(Index2d idx, float v)
{
    setV(idx.i, idx.j,v);
}

bool StaggeredVelocityGrid::getUValidity(Index2d idx)
{
    return getUValidity(idx.i, idx.j);
}

bool StaggeredVelocityGrid::getVValidity(Index2d idx)
{
    return getVValidity(idx.i, idx.j);
}

void StaggeredVelocityGrid::setUValidity(Index2d idx, bool uValidity)
{
    setUValidity(idx.i, idx.j, uValidity);
}

void StaggeredVelocityGrid::setVValidity(Index2d idx, bool vValidity)
{
    setVValidity(idx.i, idx.j, vValidity);
}

void StaggeredVelocityGrid::extrapolate(int extrapolationRadius)
{
    simmath::breadthFirstExtrapolate(m_velocityGridU,m_uSampleValidity,
                                     extrapolationRadius,
                                     m_extrapolationNeighborRadius,
                                     false);
    simmath::breadthFirstExtrapolate(m_velocityGridV,m_vSampleValidity,
                                     extrapolationRadius,
                                     m_extrapolationNeighborRadius,
                                     false);
}

Vec3 StaggeredVelocityGrid::velocityAt(float i, float j) const
{
    return Vec3(m_velocityGridU.interpolateAt(i,j),
                  m_velocityGridV.interpolateAt(i,j));
}

Vec3 StaggeredVelocityGrid::velocityAt(Vec3 position) const
{
    return velocityAt(position.x(),position.y());
}
