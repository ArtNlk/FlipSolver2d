#ifndef SDFGRID_H
#define SDFGRID_H

#include "geometry2d.h"
#include "grid2d.h"

class SdfGrid : public Grid2d<float>
{
public:
    SdfGrid(size_t sizeI, size_t sizeJ);

    Vertex closestSurfacePoint(float i, float j);
    Vertex closestSurfacePoint(Vertex pos);
};

#endif // SDFGRID_H
