#include "functions.h"

#include <cmath>

float math::frac(float v)
{
     return v-static_cast<long>(v);
}

int math::integr(float v)
{
    return static_cast<int>(std::floor(v));
}
