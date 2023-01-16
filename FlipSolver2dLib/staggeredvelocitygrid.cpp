#include "staggeredvelocitygrid.h"

#include "grid2d.h"
#include "mathfuncs.h"

StaggeredVelocityGrid::StaggeredVelocityGrid(int sizeI, int sizeJ) :
    LinearIndexable2d(sizeI, sizeJ),
    m_velocityGridU(sizeI + 1, sizeJ, 0.f, OOBStrategy::OOB_EXTEND),
    m_velocityGridV(sizeI, sizeJ + 1, 0.f, OOBStrategy::OOB_EXTEND),
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

float StaggeredVelocityGrid::getU(int i, int j)
{
    return m_velocityGridU.getAt(i,j);
}

float StaggeredVelocityGrid::getV(int i, int j)
{
    return m_velocityGridV.getAt(i,j);
}

Vertex StaggeredVelocityGrid::velocityAt(float i, float j)
{
    return Vertex(simmath::lerpUGrid(i, j, m_velocityGridU),simmath::lerpVGrid(i, j, m_velocityGridV));
}

Vertex StaggeredVelocityGrid::velocityAt(Vertex position)
{
    return velocityAt(position.x(),position.y());
}
