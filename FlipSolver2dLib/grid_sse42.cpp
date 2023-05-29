#include "grid_sse42.h"

#include "grid2d.h"
#include <array>
#include <smmintrin.h>

float Grid_sse42::cubicInterpF(const Grid2d<float> *grid, float i, float j)
{
    i+= grid->m_gridOffset.x();
    j+= grid->m_gridOffset.y();

    i -= 0.5;
    j -= 0.5;
    float iFactor = simmath::frac(i);
    float jFactor = simmath::frac(j);
    int iInt = simmath::integr(i);
    int jInt = simmath::integr(j);

    float isqrd = iFactor*iFactor;
    float iqbd = isqrd*iFactor;

    float jsqrd = jFactor*jFactor;
    float jqbd = jsqrd*jFactor;

    __m128 _iWeights;
    __m128 _jWeights;
    __m128 _half, _third, _sixth;
    __m128 _temp;
    __m128 _gridVals;
    __m128 _output;

    _half = _mm_set1_ps(1.f/2.f);
    _third = _mm_set1_ps(-1.f/3.f);
    _sixth = _mm_set1_ps(1.f/6.f);

    //I weights
    _iWeights = _mm_set_ps(0.f,0.f,1.f,0.f);
    _temp = _mm_set_ps(0.f,iFactor,-isqrd,0.f);
    _iWeights = _mm_add_ps(_iWeights,_temp);

    _temp = _mm_set_ps(0.f,isqrd-iqbd,iqbd - iFactor,isqrd);
    _temp = _mm_mul_ps(_temp,_half);
    _iWeights = _mm_add_ps(_iWeights,_temp);

    _temp = _mm_set_ps(0.f,0.f,0.f,iFactor);
    _temp = _mm_mul_ps(_temp,_third);
    _iWeights = _mm_add_ps(_iWeights,_temp);

    _temp = _mm_set_ps(iqbd - iFactor,0.f,0.f,-iqbd);
    _temp = _mm_mul_ps(_temp,_sixth);
    _iWeights = _mm_add_ps(_iWeights,_temp);

    //J weights
    _jWeights = _mm_set_ps(0.f,0.f,1.f,0.f);
    _temp = _mm_set_ps(0.f,jFactor,-jsqrd,0.f);
    _jWeights = _mm_add_ps(_jWeights,_temp);

    _temp = _mm_set_ps(0.f,jsqrd-jqbd,jqbd - jFactor,jsqrd);
    _temp = _mm_mul_ps(_temp,_half);
    _jWeights = _mm_add_ps(_jWeights,_temp);

    _temp = _mm_set_ps(0.f,0.f,0.f,jFactor);
    _temp = _mm_mul_ps(_temp,_third);
    _jWeights = _mm_add_ps(_jWeights,_temp);

    _temp = _mm_set_ps(jqbd - jFactor,0.f,0.f,-jqbd);
    _temp = _mm_mul_ps(_temp,_sixth);
    _jWeights = _mm_add_ps(_jWeights,_temp);

    _temp = _mm_set1_ps(0.f);

    int jCoord = jInt - 1;
    _gridVals = _mm_set_ps(grid->getAt(i+2,jCoord),grid->getAt(i+1,jCoord),
                           grid->getAt(i,jCoord),grid->getAt(i-1,jCoord));
    _temp = _mm_add_ps(_temp,_mm_dp_ps(_gridVals,_iWeights,0xF1));

    jCoord++;
    _gridVals = _mm_set_ps(grid->getAt(i+2,jCoord),grid->getAt(i+1,jCoord),
                           grid->getAt(i,jCoord),grid->getAt(i-1,jCoord));
    _temp = _mm_add_ps(_temp,_mm_dp_ps(_gridVals,_iWeights,0xF2));

    jCoord++;
    _gridVals = _mm_set_ps(grid->getAt(i+2,jCoord),grid->getAt(i+1,jCoord),
                           grid->getAt(i,jCoord),grid->getAt(i-1,jCoord));
    _temp = _mm_add_ps(_temp,_mm_dp_ps(_gridVals,_iWeights,0xF4));

    jCoord++;
    _gridVals = _mm_set_ps(grid->getAt(i+2,jCoord),grid->getAt(i+1,jCoord),
                           grid->getAt(i,jCoord),grid->getAt(i-1,jCoord));
    _temp = _mm_add_ps(_temp,_mm_dp_ps(_gridVals,_iWeights,0xF8));

    _output = _mm_dp_ps(_jWeights,_temp, 0xF1);

    return _mm_cvtss_f32(_output);
}
