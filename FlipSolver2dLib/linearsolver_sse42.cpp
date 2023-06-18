#include "linearsolver_sse42.h"
#include "materialgrid.h"

#include <emmintrin.h>
#include <nmmintrin.h>
#include <popcntintrin.h>
#include <smmintrin.h>

void LinearSolver_sse42::premaskPressuresThread(const Range range, const MaterialGrid& materials, std::vector<double>& pressures)
{
    const unsigned int regSize = 2;
    const unsigned int leftover = range.size() % regSize;
    const unsigned int wholeElements = range.size() - leftover;
    const unsigned int newEnd = range.start + wholeElements;
    const unsigned int leftoverEnd = newEnd + leftover;

    __m128i _solid = _mm_set1_epi64x(FluidMaterial::SOLID);

    for(unsigned int i = range.start; i < newEnd; i+=regSize)
    {
        __m128d _pressures = _mm_loadu_pd(pressures.data() + i);
        __m128i _materials = _mm_set_epi64x(materials.data()[i+1], materials.data()[i]);
        _materials = _mm_cmpeq_epi64(_materials,_solid);
        _materials = _mm_xor_si128(_materials, _mm_set1_epi64x(-1));
        _pressures = _mm_and_pd(_pressures, _mm_castsi128_pd(_materials));
        _mm_storeu_pd(pressures.data() + i,_pressures);
    }

    for(unsigned int i = newEnd; i < leftoverEnd; i++)
    {
        pressures[i] *= !solidTest(materials.data()[i]);
    }
}
