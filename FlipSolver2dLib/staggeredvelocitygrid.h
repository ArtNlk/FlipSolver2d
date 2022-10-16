#ifndef STAGGEREDVELOCITYGRID_H
#define STAGGEREDVELOCITYGRID_H

#include "linearindexable2d.h"
#include "grid2d.h"
#include "geometry2d.h"

class StaggeredVelocityGrid : public LinearIndexable2d
{
public:
    StaggeredVelocityGrid(int sizeI, int sizeJ);

    Grid2d<float>& velocityGridU();
    Grid2d<float>& velocityGridV();

    Grid2d<bool>& uSampleValidityGrid();
    Grid2d<bool>& vSampleValidityGrid();

    Vertex velocityAt(float i, float j);
    Vertex velocityAt(Vertex position);

protected:
    Grid2d<float> m_velocityGridU;
    Grid2d<float> m_velocityGridV;
    Grid2d<bool> m_uSampleValidity;
    Grid2d<bool> m_vSampleValidity;
};

#endif // STAGGEREDVELOCITYGRID_H
