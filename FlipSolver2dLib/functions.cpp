#include "functions.h"

#include <cmath>

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
