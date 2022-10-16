#include "staggeredvelocitygrid.h"

#include "mathfuncs.h"

StaggeredVelocityGrid::StaggeredVelocityGrid(int sizeI, int sizeJ) :
    LinearIndexable2d(sizeI, sizeJ),
    m_velocityGridU(sizeI + 1, sizeJ, 0.f),
    m_velocityGridV(sizeI, sizeJ + 1, 0.f),
    m_uSampleValidity(sizeI + 1, sizeJ, false),
    m_vSampleValidity(sizeI, sizeJ + 1, false)
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

Vertex StaggeredVelocityGrid::velocityAt(float i, float j)
{
    return Vertex(math::lerpUGrid(i, j, m_velocityGridU),math::lerpVGrid(i, j, m_velocityGridV));
}

Vertex StaggeredVelocityGrid::velocityAt(Vertex position)
{
    return velocityAt(position.x(),position.y());
}
