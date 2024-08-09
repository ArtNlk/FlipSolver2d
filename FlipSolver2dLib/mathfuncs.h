#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include "geometry2d.h"

template<class T>
class Grid2d;

namespace simmath
{
    float frac(float v);
    int integr(float v);
    float lerp(float a, float b, float f);
    float avg(float a,float b);
    inline float bSpline(float value)
    {
        value = std::abs(value);
        return (0.75f - value * value) * (value < 0.5f)
               + (0.5f * (1.5f - value) * (1.5f - value)) * (value >= 0.5f && value < 1.5f);
    }
    inline float quadraticBSpline(float x, float y)
    {
        return simmath::bSpline(x) * simmath::bSpline(y) * simmath::bSpline(0.f);
        //return simmath::linearHat(x) * simmath::linearHat(y);
    }
    float linearHat(float value);
    float bilinearHat(float x, float y);
    float lerpUGrid(float i, float j, const Grid2d<float> &gridU);
    float lerpVGrid(float i, float j, const Grid2d<float> &gridV);

    float lerpCenteredGrid(Vertex &position, const Grid2d<float> &grid, Vertex gridOffset = Vertex(0.f,0.f,0.f));

    float lerpCenteredGrid(float i, float j, const Grid2d<float> &grid, Vertex gridOffset = Vertex(0.f,0.f,0.f));

    void breadthFirstExtrapolate(Grid2d<float> &extrapolatedGrid, Grid2d<bool> &flagGrid,
                                          int extrapRadius, int neighborRadius, bool vonNeumannNeighborMode);

    Vertex gradCenteredGrid(ssize_t i, ssize_t j, const Grid2d<float> &grid);
}

#endif // FUNCTIONS_H
