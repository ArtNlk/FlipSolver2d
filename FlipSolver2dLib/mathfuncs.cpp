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
#include "simsettings.h"


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

float simmath::lerpUGrid(float i, float j, const Grid2d<float> &gridU)
{
    i = std::clamp(i,0.f,static_cast<float>(gridU.sizeI() - 1));
    j = std::clamp(j,0.f,static_cast<float>(gridU.sizeJ() - 1));
    Index2d currentCell(simmath::integr(i),simmath::integr(j));

    Index2d cell2(currentCell.m_i,simmath::frac(j) >= 0.5f ? currentCell.m_j + 1 : currentCell.m_j - 1);

    Index2d cell3(cell2.m_i + 1, cell2.m_j);

    Index2d cell4(currentCell.m_i + 1, currentCell.m_j);

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

    Index2d cell2(simmath::frac(i) >= 0.5f ? currentCell.m_i + 1 : currentCell.m_i - 1,currentCell.m_j);

    Index2d cell3(cell2.m_i, cell2.m_j + 1);

    Index2d cell4(currentCell.m_i, currentCell.m_j + 1);

    float iLerpFactor = simmath::frac(i) < 0.5f ? 0.5f - simmath::frac(i) : simmath::frac(i) - 0.5f;
    float jLerpFactor = simmath::frac(j);

    float v1 = simmath::lerp(gridV.getAt(currentCell),gridV.getAt(cell4),jLerpFactor);
    float v2 = simmath::lerp(gridV.getAt(cell2),gridV.getAt(cell3),jLerpFactor);

    return simmath::lerp(v1,v2,iLerpFactor);
}

Vertex simmath::gradCenteredGrid(int i, int j, const Grid2d<float> &grid)
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

void simmath::fastSweep(Grid2d<float> &values, Grid2d<bool> &extrapFlags, std::function<float (Grid2d<float> &, Vertex &, void *)> &updateFunc, void *additionalParameters)
{
    int maxIter = 10;
    float err = 0.0;
    float maxError = 1e-4f;
    for(int i = 0; i < values.sizeI(); i++)
    {
        for(int j = 0; j < values.sizeJ(); j++)
        {
            if(extrapFlags.at(i,j))
            {
                values.at(i,j) = 1e4f;
            }
        }
    }

    for(int iter = 0; iter < maxIter; iter++)
    {
        err = 0.f;
        for(int i = 0; i < values.sizeI(); i++)
        {
            for(int j = 0; j < values.sizeJ(); j++)
            {
                if(extrapFlags.at(i,j))
                {
                    Vertex pos(i,j);
                    float newValue = updateFunc(values, pos, additionalParameters);
                    if(!std::isnan(newValue))
                    {
                        err = std::max(err,std::abs(newValue-values.at(i,j)));
                        values.at(i,j) = newValue;
                    }
                }
            }
        }

        for(int i = 0; i < values.sizeI(); i++)
        {
            for(int j = values.sizeJ() - 1; j >= 0; j--)
            {
                if(extrapFlags.at(i,j))
                {
                    Vertex pos(i,j);
                    float newValue = updateFunc(values, pos,additionalParameters);
                    if(!std::isnan(newValue))
                    {
                        err = std::max(err,std::abs(newValue-values.at(i,j)));
                        values.at(i,j) = newValue;
                    }
                }
            }
        }

        for(int i = values.sizeI() - 1; i >= 0; i--)
        {
            for(int j = 0; j < values.sizeJ(); j++)
            {
                if(extrapFlags.at(i,j))
                {
                    Vertex pos(i,j);
                    float newValue = updateFunc(values, pos,additionalParameters);
                    if(!std::isnan(newValue))
                    {
                        err = std::max(err,std::abs(newValue-values.at(i,j)));
                        values.at(i,j) = newValue;
                    }
                }
            }
        }

        for(int i = values.sizeI() - 1; i >= 0; i--)
        {
            for(int j = values.sizeJ() - 1; j >= 0; j--)
            {
                if(extrapFlags.at(i,j))
                {
                    Vertex pos(i,j);
                    float newValue = updateFunc(values, pos,additionalParameters);
                    if(!std::isnan(newValue))
                    {
                        err = std::max(err,std::abs(newValue-values.at(i,j)));
                        values.at(i,j) = newValue;
                    }
                }
            }
        }

        if(err <= maxError)
        {
            std::cout << "fast sweep finished err = " << err << " iter = " << iter << '\n';
            return;
        }
    }
    std::cout << "fast sweep out of iter! err = " << err << '\n';
}

float simmath::normalDerivLinearExapolationUpdate(Grid2d<float> &grid, Vertex &pos, void* /*unused*/)
{
    int i = static_cast<int>(pos.x());
    int j = static_cast<int>(pos.y());
    float h = grid.sizeI() * grid.sizeJ();
    Vertex normal = simmath::gradCenteredGrid(i,j,grid).normalized();
    bool normalXPositive = normal.x() > 0;
    bool normalYPositive = normal.y() > 0;
    Index2d iNeighborIndex = normalXPositive? Index2d(i-1,j) : Index2d(i+1,j);
    Index2d jNeighborIndex = normalYPositive? Index2d(i,j-1) : Index2d(i,j+1);
    float aySign = normalXPositive != normalYPositive? -1 : 1;

    float ax = normal.x() * grid.getAt(iNeighborIndex);
    float ay = normal.y() * grid.getAt(jNeighborIndex) * aySign;
    float denom = normal.x() + (normal.y() * aySign);
    return (ax+ay)/denom;
}

float simmath::sdfLinearExapolationUpdate(Grid2d<float> &grid, Vertex &pos, void *normalDerivGridPtr)
{
    Grid2d<float>* normalDerivGrid = static_cast<Grid2d<float>*>(normalDerivGridPtr);
    int i = static_cast<int>(pos.x());
    int j = static_cast<int>(pos.y());
    Vertex normal = simmath::gradCenteredGrid(i,j,grid).normalized();
    bool normalXPositive = normal.x() > 0;
    bool normalYPositive = normal.y() > 0;

    Index2d iNeighborIndex = normalXPositive? Index2d(i-1,j) : Index2d(i+1,j);

    Index2d jNeighborIndex = normalYPositive? Index2d(i,j-1) : Index2d(i,j+1);

    float aySign = normalXPositive != normalYPositive? -1 : 1;
    float fSign = normalXPositive? 1 : -1;

    float ax = normal.x() * grid.getAt(iNeighborIndex);
    float ay = normal.y() * grid.getAt(jNeighborIndex) * aySign;
    float fFactor = normalDerivGrid->getAt(i,j) * fSign;
    float denom = normal.x() + (normal.y() * aySign);
    return (ax+ay+fFactor)/denom;
}

Grid2d<float> simmath::calculateCenteredGridCurvature(Grid2d<float> &grid)
{
    Grid2d<float> output(grid.sizeI(), grid.sizeJ(), 0.f, grid.oobStrat());
    Grid2d<float> tempI(grid.sizeI(), grid.sizeJ(), 0.f, OOBStrategy::OOB_EXTEND);
    Grid2d<float> tempJ(grid.sizeI(), grid.sizeJ(), 0.f, OOBStrategy::OOB_EXTEND);

//    for(int i = 0; i < grid.sizeI(); i++)
//    {
//        for(int j = 0; j < grid.sizeJ(); j++)
//        {
//            Vertex deriv = simmath::gradCenteredGrid(i,j,grid);
//            float derivI = deriv.x();
//            float derivJ = deriv.y();
//            Vertex secondDeriv = simmath::secondPartialDerivOnedir(i,j,grid);
//            float secondDerivI = secondDeriv.x();
//            float secondDerivJ = secondDeriv.y();
//            float secondDerivIj = simmath::secondPartialDerivIj(i,j,grid);

//            float firstTerm = derivI*derivI*secondDerivJ;
//            float secondTerm = 2.f*derivI*derivJ*secondDerivIj;
//            float thirdTerm = derivJ*derivJ*secondDerivI;
//            float denom = simmath::gradCenteredGrid(i,j,grid).distFromZero();
//            denom = denom*denom*denom;

//            output.at(i,j) = (firstTerm + secondTerm + thirdTerm) / denom;
//        }
//    }

    for(int i = 0; i < grid.sizeI(); i++)
    {
        for(int j = 0; j < grid.sizeJ(); j++)
        {
            Vertex grad = simmath::gradCenteredGrid(i,j,grid);
            Vertex normal = grad / grad.distFromZero();
            tempI.at(i,j) = normal.x();
            tempJ.at(i,j) = normal.y();
        }
    }

    for(int i = 0; i < grid.sizeI(); i++)
    {
        for(int j = 0; j < grid.sizeJ(); j++)
        {
            output.at(i,j) = simmath::gradCenteredGrid(i,j,tempI).x();
            output.at(i,j) += simmath::gradCenteredGrid(i,j,tempJ).y();
        }
    }

    return output;
}

Vertex simmath::secondPartialDerivOnedir(int i, int j, Grid2d<float> &grid)
{
    float currentValue = grid.at(i,j);
    float ip1Value = grid.getAt(i+1,j);
    float im1Value = grid.getAt(i-1,j);
    float jp1Value = grid.getAt(i,j+1);
    float jm1Value = grid.getAt(i,j-1);

    float derivI = ip1Value - 2*currentValue + im1Value;
    float derivJ = jp1Value - 2*currentValue + jm1Value;
    float norm = SimSettings::dx() * SimSettings::dx();

    return Vertex(derivI/norm, derivJ/norm);
}

float simmath::secondPartialDerivIj(int i, int j, Grid2d<float> &grid)
{
    float ip1jp1Value = grid.getAt(i+1,j+1);
    float im1jp1Value = grid.getAt(i-1,j+1);
    float ip1jm1Value = grid.getAt(i+1,j-1);
    float im1jm1Value = grid.getAt(i-1,j-1);

    return 0.25f * (ip1jp1Value - im1jp1Value - ip1jm1Value + im1jm1Value) / SimSettings::dx() * SimSettings::dx();
}


float simmath::lerpCenteredGrid(Vertex &position, Grid2d<float> &grid)
{
    return lerpCenteredGrid(position.x(), position.y(), grid);
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

    float v1 = simmath::lerp(grid.getAt(currentCell),grid.getAt(cell4),iLerpFactor);
    float v2 = simmath::lerp(grid.getAt(cell2),grid.getAt(cell3),iLerpFactor);

    return simmath::lerp(v1,v2,jLerpFactor);
}

void simmath::breadthFirstExtrapolate(Grid2d<float> &extrapolatedGrid, Grid2d<bool> &flagGrid,
                                      int extrapRadius, int neighborRadius, bool vonNeumannNeighborMode)
{
    int sizeI = extrapolatedGrid.sizeI();
    int sizeJ = extrapolatedGrid.sizeJ();
    Grid2d<int> markers(sizeI,sizeJ,std::numeric_limits<int>().max());
    std::queue<Index2d> wavefront;
    //Extrapolate U
    for(int i = 0; i < sizeI; i++)
    {
        for(int j = 0; j < sizeJ; j++)
        {
            if(flagGrid.at(i,j))
            {
                markers.at(i,j) = 0;
            }
        }
    }

    for(int i = 0; i < sizeI; i++)
    {
        for(int j = 0; j < sizeJ; j++)
        {
            if(markers.at(i,j) != 0)
            {
                for(Index2d& neighborIndex : extrapolatedGrid
                                                .getNeighborhood(i,j,neighborRadius,vonNeumannNeighborMode))
                {
                    if(markers.at(neighborIndex) == 0)
                    {
                        markers.at(i,j) = 1;
                        wavefront.push(Index2d(i,j));
                        break;
                    }
                }
            }
        }
    }

    while(!wavefront.empty())
    {
        Index2d index = wavefront.front();
        std::vector<Index2d> neighbors = extrapolatedGrid.getNeighborhood(index,
                                                                          neighborRadius,
                                                                          vonNeumannNeighborMode);
        double avg = 0;
        int count = 0;
        for(Index2d& neighborIndex : neighbors)
        {
            if(markers.at(neighborIndex) < markers.at(index))
            {
                avg += extrapolatedGrid.at(neighborIndex);
                count++;
            }
            if(markers.at(neighborIndex) == std::numeric_limits<int>().max() && markers.at(index) <= extrapRadius)
            {
                markers.at(neighborIndex) = markers.at(index) + 1;
                wavefront.push(neighborIndex);
            }
        }
        extrapolatedGrid.at(index) = avg / count;
        flagGrid.at(index) = true;

        wavefront.pop();
    }
}

float simmath::cubicIterpUGrid(float i, float j, const Grid2d<float> &grid)
{
    return simmath::cubicIterpGrid(i,j,grid,Vertex(0.f,-0.5f));
}

float simmath::cubicIterpVGrid(float i, float j, const Grid2d<float> &grid)
{
    return simmath::cubicIterpGrid(i,j,grid,Vertex(-0.5f,0.f));
}

float simmath::cubicIterpGrid(float i, float j, const Grid2d<float> &grid, Vertex gridOffset)
{
    i+= gridOffset.x();
    j+= gridOffset.y();
    float iFactor = simmath::frac(i);
    float jFactor = simmath::frac(j);
    float iInt = simmath::integr(i);
    float jInt = simmath::integr(j);

    auto calcWeights = [](float f)
    {
        float fSqrd = f*f;
        float fCubed = fSqrd*f;
        std::array<float,4> out = {0.f};
        out[0] = -(1.f/3.f)*f + (1.f/2.f)*fSqrd - (1.f/6.f)*fCubed;
        out[1] = 1.f - fSqrd + (1.f/2.f)*(fCubed-f);
        out[2] = f + (1.f/2.f)*(fSqrd - fCubed);
        out[3] = (1.f/6.f)*(fCubed - f);
        return out;
    };

    std::array<float,4> iWeights = calcWeights(iFactor);
    std::array<float,4> jWeights = calcWeights(jFactor);
    std::array<float,4> temp = {0.f};

    for(int index = 0; index < 4; index++)
    {
        int j = jInt + index - 1;
        temp[index] = iWeights[0]*grid.getAt(i-1,j)
                    + iWeights[1]*grid.getAt(i,j)
                    + iWeights[2]*grid.getAt(i+1,j)
                    + iWeights[3]*grid.getAt(i+2,j);
    }

    float output = jWeights[0]*temp[0]
                 + jWeights[1]*temp[1]
                 + jWeights[2]*temp[2]
                 + jWeights[3]*temp[3];
    return output;
}
