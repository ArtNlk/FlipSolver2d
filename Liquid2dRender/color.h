#ifndef COLOR_H
#define COLOR_H

class Color
{
public:
    Color() :
        m_r(0),
        m_g(0),
        m_b(0),
        m_a(0)
    {
    }

    Color(int r, int g, int b, int a = 255) :
        m_r(r),
        m_g(g),
        m_b(b),
        m_a(a)
    {
    }

    Color(float r, float g, float b, float a = 1) :
        m_r(r*255),
        m_g(g*255),
        m_b(b*255),
        m_a(a*255)
    {
    }

    inline int &r() {return m_r;}
    inline int &g() {return m_g;}
    inline int &b() {return m_b;}
    inline int &a() {return m_a;}

    inline float rf() const {return m_r / 255.f;}
    inline float gf() const {return m_g / 255.f;}
    inline float bf() const {return m_b / 255.f;}
    inline float af() const {return m_a / 255.f;}

protected:
    int m_r;
    int m_g;
    int m_b;
    int m_a;
};

#endif // COLOR_H
