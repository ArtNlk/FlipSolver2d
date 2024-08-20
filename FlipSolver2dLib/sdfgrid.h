#ifndef SDFGRID_H
#define SDFGRID_H

#include "geometry2d.h"
#include "grid2d.h"

class SdfGrid : public Grid2d<float>
{
public:
    SdfGrid(size_t sizeI, size_t sizeJ);

    Vec3 closestSurfacePoint(float i, float j);
    Vec3 closestSurfacePoint(Vec3 pos);
};

#endif // SDFGRID_H
