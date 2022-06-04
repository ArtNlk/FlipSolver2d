#include "mathfuncs.h"

#include <cmath>
#include <algorithm>

#include "customassert.h"


float math::frac(float v)
{
     return v-static_cast<long>(v);
}

int math::integr(float v)
{
    return static_cast<int>(std::floor(v));
}

float math::lerp(float a, float b, float f)
{
    return (a * (1.0f - f)) + (b * f);
}

float math::bSpline(float value)
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

float math::qudraticBSpline(float x, float y)
{
    return math::bSpline(x) * math::bSpline(y) * math::bSpline(0.f);
}

float math::lerpUGrid(float i, float j, Grid2d<float> &gridU)
{
    i = std::clamp(i,0.f,static_cast<float>(gridU.sizeI() - 1));
    j = std::clamp(j,0.f,static_cast<float>(gridU.sizeJ() - 1));
    Index2d currentCell(math::integr(i),math::integr(j));

    Index2d cell2(currentCell.m_i,math::frac(j) >= 0.5f ? currentCell.m_j + 1 : currentCell.m_j - 1);

    Index2d cell3(cell2.m_i + 1, cell2.m_j);

    Index2d cell4(currentCell.m_i + 1, currentCell.m_j);

    float iLerpFactor = math::frac(i);
    float jLerpFactor = math::frac(j) < 0.5f ? 0.5f - math::frac(j) : math::frac(j) - 0.5f;

    float v1 = gridU.inBounds(cell4) ? math::lerp(gridU.at(currentCell),gridU.at(cell4),iLerpFactor) : gridU.at(currentCell);
    float v2 = gridU.inBounds(cell3) ? math::lerp(gridU.at(cell2),gridU.at(cell3),iLerpFactor) : v1;

    return math::lerp(v1,v2,jLerpFactor);
}

float math::lerpVGrid(float i, float j, Grid2d<float> &gridV)
{
    i = std::clamp(i,0.f,static_cast<float>(gridV.sizeI() - 1));
    j = std::clamp(j,0.f,static_cast<float>(gridV.sizeJ() - 1));
    Index2d currentCell(math::integr(i),math::integr(j));

    Index2d cell2(math::frac(i) >= 0.5f ? currentCell.m_i + 1 : currentCell.m_i - 1,currentCell.m_j);

    Index2d cell3(cell2.m_i, cell2.m_j + 1);

    Index2d cell4(currentCell.m_i, currentCell.m_j + 1);

    float iLerpFactor = math::frac(i) < 0.5f ? 0.5f - math::frac(i) : math::frac(i) - 0.5f;
    float jLerpFactor = math::frac(j);

    float v1 = gridV.inBounds(cell4) ? math::lerp(gridV.at(currentCell),gridV.at(cell4),jLerpFactor) : gridV.at(currentCell);
    float v2 = gridV.inBounds(cell3) ? math::lerp(gridV.at(cell2),gridV.at(cell3),jLerpFactor) : v1;

    return math::lerp(v1,v2,iLerpFactor);
}

float math::lerpCenteredGrid(float i, float j, Grid2d<float> &grid)
{
    i = std::clamp(i,0.f,static_cast<float>(grid.sizeI() - 1));
    j = std::clamp(j,0.f,static_cast<float>(grid.sizeJ() - 1));
    Index2d currentCell(math::integr(i),math::integr(j));

    Index2d cell2(currentCell.m_i,math::frac(j) >= 0.5f ? currentCell.m_j + 1 : currentCell.m_j - 1);

    Index2d cell3(math::frac(i) >= 0.5f ? currentCell.m_i + 1 : currentCell.m_i - 1, math::frac(j) >= 0.5f ? currentCell.m_j + 1 : currentCell.m_j - 1);

    Index2d cell4(math::frac(i) >= 0.5f ? currentCell.m_i + 1 : currentCell.m_i - 1, currentCell.m_j);

    float iLerpFactor = math::frac(i) < 0.5f ? 0.5f - math::frac(i) : math::frac(i) - 0.5f;
    float jLerpFactor = math::frac(j) < 0.5f ? 0.5f - math::frac(j) : math::frac(j) - 0.5f;

    float v1 = grid.inBounds(cell4) ? math::lerp(grid.at(currentCell),grid.at(cell4),iLerpFactor) : grid.at(currentCell);
    float v2 = grid.inBounds(cell3) ? math::lerp(grid.at(cell2),grid.at(cell3),iLerpFactor) : v1;

    return math::lerp(v1,v2,jLerpFactor);
}

Vertex math::gradCenteredGrid(float i, float j, Grid2d<float> &grid)
{
    Index2d currentCell(math::integr(i),math::integr(j));
    if(!grid.inBounds(currentCell)) return 0.f;

    Index2d neighborI(math::frac(i) >= 0.5f ? currentCell.m_i + 1 : currentCell.m_i - 1, currentCell.m_j);
    if(!grid.inBounds(neighborI)) neighborI.m_i = currentCell.m_i + -1*(neighborI.m_i - currentCell.m_i);

    Index2d neighborJ(currentCell.m_i,math::frac(j) >= 0.5f ? currentCell.m_j + 1 : currentCell.m_j - 1);
    if(!grid.inBounds(neighborJ)) neighborJ.m_j = currentCell.m_j + -1*(neighborJ.m_j - currentCell.m_j);

    float gradI = math::frac(j) >= 0.5f ? grid.at(currentCell) - grid.at(neighborI) : grid.at(neighborI) - grid.at(currentCell);
    float gradJ = math::frac(i) >= 0.5f ? grid.at(currentCell) - grid.at(neighborJ) : grid.at(neighborJ) - grid.at(currentCell);

    return Vertex(gradI, gradJ);
}
