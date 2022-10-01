#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include "grid2d.h"
#include "geometry2d.h"

namespace math
{
    float frac(float v);
    int integr(float v);
    float lerp(float a, float b, float f);
    float avg(float a,float b);
    float bSpline(float value);
    float quadraticBSpline(float x, float y);
    float linearHat(float value);
    float bilinearHat(float x, float y);
    float lerpUGrid(float i, float j, Grid2d<float> &gridU);
    float lerpVGrid(float i, float j, Grid2d<float> &gridV);
    float lerpCenteredGrid(float i, float j, Grid2d<float> &grid);
    Vertex gradCenteredGrid(float i, float j, Grid2d<float> &grid);
}

#endif // FUNCTIONS_H
