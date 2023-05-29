#ifndef GRID_SSE42_H
#define GRID_SSE42_H

#include <nmmintrin.h>

#include "mathfuncs.h"

template<class T>
class Grid2d;

class Grid_sse42
{
public:
    Grid_sse42() = default;

    static float cubicInterpF(const Grid2d<float>* grid, float i, float j);
};

#endif // GRID_SSE42_H
