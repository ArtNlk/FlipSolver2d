#include "linearsolver_sse42.h"

#include <emmintrin.h>
#include <nmmintrin.h>
#include <popcntintrin.h>
#include <smmintrin.h>

#include "linearsolver.h"

void LinearSolver_sse42::dampedJacobiThread(LinearSolver*, const Range range, const MaterialGrid &materials, std::vector<double> &vout, const std::vector<double> &pressures, const std::vector<double> &rhs)
{
    const double tune = 2.0/3.0;
    for(int i = range.start; i < range.end; i++)
    {
        if(!fluidTest(materials.data()[i]))
        {
            continue;
        }
        const auto weights = getMultigridMatrixEntriesForCell(materials,i);
        __m128d _weightsI, _weightsJ;
        __m128d _pressuresI, _pressuresJ;
        __m128d _temp;
        _weightsI = _mm_loadu_pd(&weights.second[0]);
        _weightsJ = _mm_loadu_pd(&weights.second[2]);
        _pressuresI = _mm_setr_pd(pressures[weights.first[0]],
                                 pressures[weights.first[1]]);
        _pressuresJ = _mm_setr_pd(pressures[weights.first[2]],
                                 pressures[weights.first[3]]);
        _temp = _mm_dp_pd(_pressuresI,_weightsI,0xF1);
        _temp = _mm_add_pd(_temp,_mm_dp_pd(_pressuresJ,_weightsJ,0xF2));
        double result = _temp[0] + _temp[1];
        vout[i] = ((rhs[i]-result)/-4.0) * tune;
    }
}

std::pair<std::array<int,4>,std::array<double,4>> LinearSolver_sse42::getMultigridMatrixEntriesForCell(const MaterialGrid &materials, int linearIdx)
{
    std::pair<std::array<int,4>,std::array<double,4>> output;
    std::array<int,4> neighbors = materials.immidiateNeighbors(linearIdx);
    const std::vector<FluidMaterial> &materialsData = materials.data();
    const auto dataSize = materialsData.size();

    int outputIdx = 0;
    for(int idx : neighbors)
    {
        if(idx < 0 || idx >= dataSize)
        {
            output.first[outputIdx] = 0;
            output.second[outputIdx] = 0.;
            outputIdx++;
            continue;
        }
        const double weight = solidTest(materialsData[idx]) ? 0. : 1.;
        output.first[outputIdx] = idx;
        output.second[outputIdx] = weight;
        outputIdx++;
    }

    return output;
}
