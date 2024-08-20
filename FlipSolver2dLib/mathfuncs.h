#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <functional>
#include <queue>
#include <type_traits>

#include "geometry2d.h"
#include "index2d.h"

template<class T>
class Grid2d;

namespace simmath
{
    float frac(float v);
    int integr(float v);
    float lerp(float a, float b, float f);
    float avg(float a,float b);
    float bSpline(float value);
    float quadraticBSpline(float x, float y);
    float linearHat(float value);
    float bilinearHat(float x, float y);
    float lerpUGrid(float i, float j, const Grid2d<float> &gridU);
    float lerpVGrid(float i, float j, const Grid2d<float> &gridV);

    float lerpCenteredGrid(Vec3 &position, const Grid2d<float> &grid, Vec3 gridOffset = Vec3(0.f,0.f,0.f));

    float lerpCenteredGrid(float i, float j, const Grid2d<float> &grid, Vec3 gridOffset = Vec3(0.f,0.f,0.f));

    void breadthFirstExtrapolate(Grid2d<float> &extrapolatedGrid, Grid2d<bool> &flagGrid,
                                          int extrapRadius, int neighborRadius, bool vonNeumannNeighborMode);

    Vec3 gradCenteredGrid(ssize_t i, ssize_t j, const Grid2d<float> &grid);
}

#endif // FUNCTIONS_H
