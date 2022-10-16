#include "color.h"
#include "mathfuncs.h"

Color Color::fromHSVA(float h, float s, float v, float a)
{
    float r = 0, g = 0, b = 0;

    if (s == 0)
    {
        r = v;
        g = v;
        b = v;
    }
    else
    {
        int i;
        float f, p, q, t;

        if (h == 1.f)
            h = 0.f;
        else
            h = h / 0.166f;

        i = static_cast<int>(trunc(h));
        f = h - i;

        p = v * (1.0 - s);
        q = v * (1.0 - (s * f));
        t = v * (1.0 - (s * (1.0 - f)));

        switch (i)
        {
        case 0:
            r = v;
            g = t;
            b = p;
            break;

        case 1:
            r = q;
            g = v;
            b = p;
            break;

        case 2:
            r = p;
            g = v;
            b = t;
            break;

        case 3:
            r = p;
            g = q;
            b = v;
            break;

        case 4:
            r = t;
            g = p;
            b = v;
            break;

        default:
            r = v;
            g = p;
            b = q;
            break;
        }

    }

    Color c(r,g,b);

    return c;
}

Color Color::lerp(Color first, Color second, float factor)
{
    float r = math::lerp(first.rf(),second.rf(),factor);
    float g = math::lerp(first.gf(),second.gf(),factor);
    float b = math::lerp(first.bf(),second.bf(),factor);
    float a = math::lerp(first.af(),second.af(),factor);

    return Color(r,g,b,a);
}
