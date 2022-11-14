#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <functional>

#include "grid2d.h"
#include "geometry2d.h"

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
    float lerpUGrid(float i, float j, Grid2d<float> &gridU);
    float lerpVGrid(float i, float j, Grid2d<float> &gridV);
    float lerpCenteredGrid(float i, float j, Grid2d<float> &grid);
    Vertex gradCenteredGrid(int i, int j, Grid2d<float> &grid);

    void fastSweep(Grid2d<float> &values,
                   Grid2d<bool> &extrapFlags,
                   std::function<float(Grid2d<float>&,Vertex&, void*)> &updateFunc,
                   void* additionalParameters);

    float normalDerivLinearExapolationUpdate(Grid2d<float> &grid, Vertex& pos, void*);
    float sdfLinearExapolationUpdate(Grid2d<float> &grid, Vertex& pos, void* normalDerivGrid);
}

#endif // FUNCTIONS_H
