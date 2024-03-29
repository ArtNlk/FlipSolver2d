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

    float lerpCenteredGrid(Vertex &position, const Grid2d<float> &grid, Vertex gridOffset = Vertex(0.f,0.f,0.f));

    float lerpCenteredGrid(float i, float j, const Grid2d<float> &grid, Vertex gridOffset = Vertex(0.f,0.f,0.f));

    void breadthFirstExtrapolate(Grid2d<float> &extrapolatedGrid, Grid2d<bool> &flagGrid,
                                          int extrapRadius, int neighborRadius, bool vonNeumannNeighborMode);

    Vertex gradCenteredGrid(int i, int j, const Grid2d<float> &grid);

    void fastSweep(Grid2d<float> &values,
                   Grid2d<bool> &extrapFlags,
                   std::function<float(Grid2d<float>&,Vertex&, void*)> &updateFunc,
                   void* additionalParameters);

    float normalDerivLinearExapolationUpdate(Grid2d<float> &grid, Vertex& pos, void*);
    float sdfLinearExapolationUpdate(Grid2d<float> &grid, Vertex& pos, void* normalDerivGrid);
    Grid2d<float> calculateCenteredGridCurvature(Grid2d<float>& grid);
    Vertex secondPartialDerivOnedir(int i, int j, Grid2d<float> &grid, float dx);
    float secondPartialDerivIj(int i, int j, Grid2d<float> &grid, float dx);
}

#endif // FUNCTIONS_H
