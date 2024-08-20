#include "mathfuncs.h"

#include <array>
#include <cmath>
#include <algorithm>
#include <limits>
#include <queue>
#include <type_traits>

#include "customassert.h"
#include "grid2d.h"
#include "index2d.h"
#include "logger.h"



float simmath::frac(float v)
{
     return v-static_cast<long>(v);
}

int simmath::integr(float v)
{
    return static_cast<int>(std::floor(v));
}

float simmath::lerp(float a, float b, float f)
{
    return (a * (1.0f - f)) + (b * f);
}

float simmath::bSpline(float value)
{
    value = std::abs(value);
    return (0.75f - value * value)*(value < 0.5f) +
        (0.5f * (1.5f - value) * (1.5f - value))*(value >= 0.5f && value < 1.5f);
}

float simmath::quadraticBSpline(float x, float y)
{
    return simmath::bSpline(x) * simmath::bSpline(y) * simmath::bSpline(0.f);
    //return simmath::linearHat(x) * simmath::linearHat(y);
}

float simmath::lerpUGrid(float i, float j, const Grid2d<float> &gridU)
{
    i = std::clamp(i,0.f,static_cast<float>(gridU.sizeI() - 1));
    j = std::clamp(j,0.f,static_cast<float>(gridU.sizeJ() - 1));
    Index2d currentCell(simmath::integr(i),simmath::integr(j));

    Index2d cell2(currentCell.i,simmath::frac(j) >= 0.5f ? currentCell.j + 1 : currentCell.j - 1);

    Index2d cell3(cell2.i + 1, cell2.j);

    Index2d cell4(currentCell.i + 1, currentCell.j);

    float iLerpFactor = simmath::frac(i);
    float jLerpFactor = simmath::frac(j) < 0.5f ? 0.5f - simmath::frac(j) : simmath::frac(j) - 0.5f;

    float v1 = simmath::lerp(gridU.getAt(currentCell),gridU.getAt(cell4),iLerpFactor);
    float v2 = simmath::lerp(gridU.getAt(cell2),gridU.getAt(cell3),iLerpFactor);

    return simmath::lerp(v1,v2,jLerpFactor);
}

float simmath::lerpVGrid(float i, float j, const Grid2d<float> &gridV)
{
    i = std::clamp(i,0.f,static_cast<float>(gridV.sizeI() - 1));
    j = std::clamp(j,0.f,static_cast<float>(gridV.sizeJ() - 1));
    Index2d currentCell(simmath::integr(i),simmath::integr(j));

    Index2d cell2(simmath::frac(i) >= 0.5f ? currentCell.i + 1 : currentCell.i - 1,currentCell.j);

    Index2d cell3(cell2.i, cell2.j + 1);

    Index2d cell4(currentCell.i, currentCell.j + 1);

    float iLerpFactor = simmath::frac(i) < 0.5f ? 0.5f - simmath::frac(i) : simmath::frac(i) - 0.5f;
    float jLerpFactor = simmath::frac(j);

    float v1 = simmath::lerp(gridV.getAt(currentCell),gridV.getAt(cell4),jLerpFactor);
    float v2 = simmath::lerp(gridV.getAt(cell2),gridV.getAt(cell3),jLerpFactor);

    return simmath::lerp(v1,v2,iLerpFactor);
}

Vertex simmath::gradCenteredGrid(ssize_t i, ssize_t j, const Grid2d<float> &grid)
{
    float neighborIp1 = grid.getAt(i+1,j);
    float neighborIm1 = grid.getAt(i-1,j);
    float neighborJp1 = grid.getAt(i,j+1);
    float neighborJm1 = grid.getAt(i,j-1);

    float gradI = neighborIp1 - neighborIm1;
    float gradJ = neighborJp1 - neighborJm1;

    return Vertex(gradI/2.f, gradJ/2.f);
}

float simmath::linearHat(float value)
{
    if(value >= 0 && value <= 1)
    {
        return 1-value;
    }
    else if(value < 0 && value >= -1)
    {
        return 1+value;
    }
    return 0;
}

float simmath::bilinearHat(float x, float y)
{
    return simmath::linearHat(x) * simmath::linearHat(y) * simmath::linearHat(0.f);
}

float simmath::avg(float a, float b)
{
    return (a+b)/2.f;
}

float simmath::lerpCenteredGrid(Vertex &position, const Grid2d<float> &grid, Vertex gridOffset)
{
    return lerpCenteredGrid(position.x(), position.y(), grid, gridOffset);
}

float simmath::lerpCenteredGrid(float i, float j, const Grid2d<float> &grid, Vertex gridOffset)
{
    i+= gridOffset.x();
    j+= gridOffset.y();

    i = std::clamp(i,0.f,static_cast<float>(grid.sizeI() - 1));
    j = std::clamp(j,0.f,static_cast<float>(grid.sizeJ() - 1));
    Index2d currentCell(simmath::integr(i),simmath::integr(j));

    Index2d cell2(currentCell.i,simmath::frac(j) >= 0.5f ? currentCell.j + 1 : currentCell.j - 1);

    Index2d cell3(simmath::frac(i) >= 0.5f ? currentCell.i + 1 : currentCell.i - 1, simmath::frac(j) >= 0.5f ? currentCell.j + 1 : currentCell.j - 1);

    Index2d cell4(simmath::frac(i) >= 0.5f ? currentCell.i + 1 : currentCell.i - 1, currentCell.j);

    float iLerpFactor = simmath::frac(i) < 0.5f ? 0.5f - simmath::frac(i) : simmath::frac(i) - 0.5f;
    float jLerpFactor = simmath::frac(j) < 0.5f ? 0.5f - simmath::frac(j) : simmath::frac(j) - 0.5f;

    float v1 = simmath::lerp(grid.getAt(currentCell),grid.getAt(cell4),iLerpFactor);
    float v2 = simmath::lerp(grid.getAt(cell2),grid.getAt(cell3),iLerpFactor);

    return simmath::lerp(v1,v2,jLerpFactor);
}

void simmath::breadthFirstExtrapolate(Grid2d<float> &extrapolatedGrid, Grid2d<bool> &flagGrid,
                                      int extrapRadius, int neighborRadius, bool vonNeumannNeighborMode)
{
    size_t sizeI = extrapolatedGrid.sizeI();
    size_t sizeJ = extrapolatedGrid.sizeJ();
    Grid2d<int> markers(sizeI,sizeJ,std::numeric_limits<int>().max());
    std::queue<size_t> wavefront;
    for(size_t i = 0; i < sizeI; i++)
    {
        for(size_t j = 0; j < sizeJ; j++)
        {
            if(flagGrid.at(i,j))
            {
                markers.at(i,j) = 0;
            }
        }
    }

    for(size_t i = 0; i < sizeI; i++)
    {
        for(size_t j = 0; j < sizeJ; j++)
        {
            if(markers.at(i,j) != 0)
            {
                for(ssize_t neighborIndex : extrapolatedGrid.getNeighborhood(i,j))
                {
                    if(neighborIndex == -1) continue;
                    if(markers.data().at(neighborIndex) == 0)
                    {
                        markers.at(i,j) = 1;
                        wavefront.push(flagGrid.linearIndex(i,j));
                        break;
                    }
                }
            }
        }
    }

    while(!wavefront.empty())
    {
        size_t index = wavefront.front();
        std::array<ssize_t, 8> neighbors = extrapolatedGrid.getNeighborhood(index);
        double avg = 0;
        int count = 0;
        for(ssize_t neighborIndex : neighbors)
        {
            if(neighborIndex == -1) continue;
            if(markers.data().at(neighborIndex) < markers.data().at(index))
            {
                avg += extrapolatedGrid.data().at(neighborIndex);
                count++;
            }
            if(markers.data().at(neighborIndex) == std::numeric_limits<int>().max() && markers.data().at(index) <= extrapRadius)
            {
                markers.data().at(neighborIndex) = markers.data().at(index) + 1;
                wavefront.push(neighborIndex);
            }
        }
        extrapolatedGrid.data().at(index) = avg / count;
        flagGrid.data().at(index) = true;

        wavefront.pop();
    }
}
