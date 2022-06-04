#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include "grid2d.h"
#include "geometry2d.h"

namespace math
{
    float frac(float v);
    int integr(float v);
    float lerp(float a, float b, float f);
    float bSpline(float value);
    float qudraticBSpline(float x, float y);
    float lerpUGrid(float i, float j, Grid2d<float> &gridU);
    float lerpVGrid(float i, float j, Grid2d<float> &gridV);
    float lerpCenteredGrid(float i, float j, Grid2d<float> &grid);
    Vertex gradCenteredGrid(float i, float j, Grid2d<float> &grid);
}

#endif // FUNCTIONS_H
