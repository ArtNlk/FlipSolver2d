#include "staggeredvelocitygrid.h"

#include "grid2d.h"
#include "mathfuncs.h"

StaggeredVelocityGrid::StaggeredVelocityGrid(int sizeI, int sizeJ) :
    LinearIndexable2d(sizeI, sizeJ),
    m_velocityGridU(sizeI + 1, sizeJ, 0.f, OOBStrategy::OOB_EXTEND, 0.f, Vertex(0.5f,0.f)),
    m_velocityGridV(sizeI, sizeJ + 1, 0.f, OOBStrategy::OOB_EXTEND, 0.f, Vertex(0.f, 0.5f)),
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

Grid2d<bool> &StaggeredVelocityGrid::uSampleValidityGrid()
{
    return m_uSampleValidity;
}

Grid2d<bool> &StaggeredVelocityGrid::vSampleValidityGrid()
{
    return m_vSampleValidity;
}

float &StaggeredVelocityGrid::u(int i, int j)
{
    return m_velocityGridU.at(i,j);
}

float &StaggeredVelocityGrid::v(int i, int j)
{
    return m_velocityGridV.at(i,j);
}

float StaggeredVelocityGrid::getU(int i, int j) const
{
    return m_velocityGridU.getAt(i,j);
}

float StaggeredVelocityGrid::getV(int i, int j) const
{
    return m_velocityGridV.getAt(i,j);
}

void StaggeredVelocityGrid::setU(int i, int j, float u)
{
    m_velocityGridU.setAt(i,j,u);
}

void StaggeredVelocityGrid::setV(int i, int j, float v)
{
    m_velocityGridV.setAt(i,j,v);
}

bool StaggeredVelocityGrid::getUValidity(int i, int j)
{
    return m_uSampleValidity.getAt(i,j);
}

bool StaggeredVelocityGrid::getVValidity(int i, int j)
{
    return m_vSampleValidity.getAt(i,j);
}

void StaggeredVelocityGrid::setUValidity(int i, int j, bool uValidity)
{
    m_uSampleValidity.setAt(i,j, uValidity);
}

void StaggeredVelocityGrid::setVValidity(int i, int j, bool vValidity)
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

Vertex StaggeredVelocityGrid::velocityAt(float i, float j) const
{
    return Vertex(m_velocityGridU.interpolateAt(i,j),
                  m_velocityGridV.interpolateAt(i,j));
}

Vertex StaggeredVelocityGrid::velocityAt(Vertex position) const
{
    return velocityAt(position.x(),position.y());
}
