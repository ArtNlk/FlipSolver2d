#include "mathfuncs.h"

#include <cmath>
#include <algorithm>

#include "customassert.h"


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
    if(value >= -1.5f && value < -0.5f)
    {
        return 0.5f*(value + 1.5f)*(value + 1.5f);
    }
    else if(value >= -0.5 && value < 0.5f)
    {
        return 0.75f-value*value;
    }
    else if(value >= 0.5f && value < 1.5f)
    {
        return 0.5f*(1.5f - value)*(1.5f - value);
    }

    return 0;
}

float simmath::quadraticBSpline(float x, float y)
{
    return simmath::bSpline(x) * simmath::bSpline(y) * simmath::bSpline(0.f);
}

float simmath::lerpUGrid(float i, float j, Grid2d<float> &gridU)
{
    i = std::clamp(i,0.f,static_cast<float>(gridU.sizeI() - 1));
    j = std::clamp(j,0.f,static_cast<float>(gridU.sizeJ() - 1));
    Index2d currentCell(simmath::integr(i),simmath::integr(j));

    Index2d cell2(currentCell.m_i,simmath::frac(j) >= 0.5f ? currentCell.m_j + 1 : currentCell.m_j - 1);

    Index2d cell3(cell2.m_i + 1, cell2.m_j);

    Index2d cell4(currentCell.m_i + 1, currentCell.m_j);

    float iLerpFactor = simmath::frac(i);
    float jLerpFactor = simmath::frac(j) < 0.5f ? 0.5f - simmath::frac(j) : simmath::frac(j) - 0.5f;

    float v1 = gridU.inBounds(cell4) ? simmath::lerp(gridU.at(currentCell),gridU.at(cell4),iLerpFactor) : gridU.at(currentCell);
    float v2 = gridU.inBounds(cell3) ? simmath::lerp(gridU.at(cell2),gridU.at(cell3),iLerpFactor) : v1;

    return simmath::lerp(v1,v2,jLerpFactor);
}

float simmath::lerpVGrid(float i, float j, Grid2d<float> &gridV)
{
    i = std::clamp(i,0.f,static_cast<float>(gridV.sizeI() - 1));
    j = std::clamp(j,0.f,static_cast<float>(gridV.sizeJ() - 1));
    Index2d currentCell(simmath::integr(i),simmath::integr(j));

    Index2d cell2(simmath::frac(i) >= 0.5f ? currentCell.m_i + 1 : currentCell.m_i - 1,currentCell.m_j);

    Index2d cell3(cell2.m_i, cell2.m_j + 1);

    Index2d cell4(currentCell.m_i, currentCell.m_j + 1);

    float iLerpFactor = simmath::frac(i) < 0.5f ? 0.5f - simmath::frac(i) : simmath::frac(i) - 0.5f;
    float jLerpFactor = simmath::frac(j);

    float v1 = gridV.inBounds(cell4) ? simmath::lerp(gridV.at(currentCell),gridV.at(cell4),jLerpFactor) : gridV.at(currentCell);
    float v2 = gridV.inBounds(cell3) ? simmath::lerp(gridV.at(cell2),gridV.at(cell3),jLerpFactor) : v1;

    return simmath::lerp(v1,v2,iLerpFactor);
}

float simmath::lerpCenteredGrid(float i, float j, Grid2d<float> &grid)
{
    i = std::clamp(i,0.f,static_cast<float>(grid.sizeI() - 1));
    j = std::clamp(j,0.f,static_cast<float>(grid.sizeJ() - 1));
    Index2d currentCell(simmath::integr(i),simmath::integr(j));

    Index2d cell2(currentCell.m_i,simmath::frac(j) >= 0.5f ? currentCell.m_j + 1 : currentCell.m_j - 1);

    Index2d cell3(simmath::frac(i) >= 0.5f ? currentCell.m_i + 1 : currentCell.m_i - 1, simmath::frac(j) >= 0.5f ? currentCell.m_j + 1 : currentCell.m_j - 1);

    Index2d cell4(simmath::frac(i) >= 0.5f ? currentCell.m_i + 1 : currentCell.m_i - 1, currentCell.m_j);

    float iLerpFactor = simmath::frac(i) < 0.5f ? 0.5f - simmath::frac(i) : simmath::frac(i) - 0.5f;
    float jLerpFactor = simmath::frac(j) < 0.5f ? 0.5f - simmath::frac(j) : simmath::frac(j) - 0.5f;

    float v1 = grid.inBounds(cell4) ? simmath::lerp(grid.at(currentCell),grid.at(cell4),iLerpFactor) : grid.at(currentCell);
    float v2 = grid.inBounds(cell3) ? simmath::lerp(grid.at(cell2),grid.at(cell3),iLerpFactor) : v1;

    return simmath::lerp(v1,v2,jLerpFactor);
}

Vertex simmath::gradCenteredGrid(float i, float j, Grid2d<float> &grid)
{
    Index2d currentCell(simmath::integr(i),simmath::integr(j));
    if(!grid.inBounds(currentCell)) return 0.f;

    Index2d neighborI(simmath::frac(i) >= 0.5f ? currentCell.m_i + 1 : currentCell.m_i - 1, currentCell.m_j);
    if(!grid.inBounds(neighborI)) neighborI.m_i = currentCell.m_i + -1*(neighborI.m_i - currentCell.m_i);

    Index2d neighborJ(currentCell.m_i,simmath::frac(j) >= 0.5f ? currentCell.m_j + 1 : currentCell.m_j - 1);
    if(!grid.inBounds(neighborJ)) neighborJ.m_j = currentCell.m_j + -1*(neighborJ.m_j - currentCell.m_j);

    float gradI = simmath::frac(j) >= 0.5f ? grid.at(currentCell) - grid.at(neighborI) : grid.at(neighborI) - grid.at(currentCell);
    float gradJ = simmath::frac(i) >= 0.5f ? grid.at(currentCell) - grid.at(neighborJ) : grid.at(neighborJ) - grid.at(currentCell);

    return Vertex(gradI, gradJ);
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
