#include "linearsolver_sse42.h"

#include <emmintrin.h>
#include <nmmintrin.h>
#include <popcntintrin.h>
#include <smmintrin.h>

#include "linearsolver.h"

void LinearSolver_sse42::dampedJacobiThread(LinearSolver *, const Range range, const MaterialGrid &materials, std::vector<double> &vout, const std::vector<double> &pressures, const std::vector<double> &rhs)
{
    const float tune = 2.f/3.f;
    const std::vector<FluidMaterial> &materialsData = materials.data();
    const auto dataSize = materialsData.size();

    __m128i _lidxMax, _lidxMin;
    _lidxMax = _mm_set1_epi64x(dataSize);
    _lidxMin = _mm_set1_epi64x(0);

    int64_t ip1LinOfst = materials.sizeJ();
    int64_t im1LinOfst = -1*ip1LinOfst;
    int64_t jp1LinOfst = 1;
    int64_t jm1LinOfst = -1;

    for(int i = range.start; i < range.end; i++)
    {
        //const auto weights = solver->getMultigridMatrixEntriesForCell(materials,i);

        __m128i _neighborIdxsI, _neighborOffsetsI, _solidityI, _validI;
        __m128d _pressuresI;
        __m128d _pressuresJ;
        __m128i _neighborIdxsJ, _neighborOffsetsJ, _solidityJ, _validJ;
        __m128i _temp, _solidityMask;

        _pressuresI = _mm_setzero_pd();
        _pressuresJ = _mm_setzero_pd();

        _solidityMask = _mm_set1_epi64x(FluidMaterial::SOLID);

        _neighborOffsetsI = _mm_set_epi64x(ip1LinOfst,im1LinOfst);
        _neighborOffsetsJ = _mm_set_epi64x(jp1LinOfst,jm1LinOfst);

        _neighborIdxsI = _mm_set1_epi64x(i);
        _neighborIdxsJ = _mm_set1_epi64x(i);

        _neighborIdxsI = _mm_add_epi64(_neighborIdxsI,_neighborOffsetsI);
        _neighborIdxsJ = _mm_add_epi64(_neighborIdxsJ,_neighborOffsetsJ);

        _temp = _mm_cmpgt_epi64(_neighborIdxsI,_lidxMin);
        _validI = _mm_cmpgt_epi64(_lidxMax,_neighborIdxsI);
        _validI = _mm_and_si128(_temp, _validI);

        _temp = _mm_cmpgt_epi64(_neighborIdxsJ,_lidxMin);
        _validJ = _mm_cmpgt_epi64(_lidxMax,_neighborIdxsJ);
        _validJ = _mm_and_si128(_temp, _validJ);

        _solidityI = _mm_set1_epi64x(materials.oobVal());
        _solidityJ = _mm_set1_epi64x(materials.oobVal());

        //No gather for sse(
        if(_validI[0])
        {
            _solidityI = _mm_insert_epi64(_solidityI,materialsData[_neighborIdxsI[0]],0);
            _pressuresI = _mm_set1_pd(pressures[_neighborIdxsI[0]]);
        }
        if(_validI[1])
        {
            _solidityI = _mm_insert_epi64(_solidityI,materialsData[_neighborIdxsI[1]],1);
            _pressuresI = _mm_blend_pd(_pressuresI,_mm_set1_pd(pressures[_neighborIdxsI[1]]),1);
        }
        if(_validJ[0])
        {
            _solidityJ = _mm_insert_epi64(_solidityJ,materialsData[_neighborIdxsJ[0]],0);
            _pressuresJ = _mm_set1_pd(pressures[_neighborIdxsJ[0]]);
        }
        if(_validJ[1])
        {
            _solidityJ = _mm_insert_epi64(_solidityJ,materialsData[_neighborIdxsJ[1]],1);
            _pressuresJ = _mm_blend_pd(_pressuresJ,_mm_set1_pd(pressures[_neighborIdxsJ[1]]),1);
        }

        //_neighborSolidity = _mm_and_si128(_neighborSolidity,_solidityMask);
        _solidityI = _mm_cmpeq_epi8(_solidityI,_solidityMask);
        _solidityJ = _mm_cmpeq_epi8(_solidityJ,_solidityMask);
        uint16_t solids = _mm_movemask_pd(_mm_castsi128_pd(_solidityI));
        solids <<= 2;
        solids |= _mm_movemask_pd(_mm_castsi128_pd(_solidityJ));
        double solidNeighborCount = static_cast<double>(_mm_popcnt_u32(solids));

        //        int currIdx = weights[0].first;
        //        float result = 0.f;
        //        for(int wIdx = 1; wIdx < weights.size(); wIdx++)
        //        {
        //            if(weights[wIdx].first < 0)
        //            {
        //                continue;
        //            }

        //            result += weights[wIdx].second * pressures[weights[wIdx].first];
        //        }
        //        vout[currIdx] = ((rhs[currIdx]-result)/weights[0].second) * tune;

        _pressuresI = _mm_castsi128_pd(_mm_and_si128(_mm_castpd_si128(_pressuresI),_solidityI));
        _pressuresJ = _mm_castsi128_pd(_mm_and_si128(_mm_castpd_si128(_pressuresJ),_solidityJ));

        _pressuresI = _mm_hadd_pd(_pressuresI, _pressuresJ); //reusing PressureI as pd temp

        double result = _pressuresI[0] + _pressuresI[1];

        vout[i] = ((rhs[i] - result)/(-solidNeighborCount)) * tune;
    }
}
