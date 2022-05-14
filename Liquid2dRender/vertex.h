#ifndef VERTEX_H
#define VERTEX_H

#include <cmath>

class Vertex
{
public:
    Vertex(float x = 0.0f, float y = 0.0f, float z = 0.0f) :
        m_x(x),
        m_y(y),
        m_z(z)
    {
    }

    inline float &x() {return m_x;}
    inline float &y() {return m_y;}
    inline float &z() {return m_z;}

    inline float distFromZero() {return std::sqrt(m_x * m_x + m_y * m_y + m_z * m_z);}

protected:
    float m_x;
    float m_y;
    float m_z;
};

#endif // VERTEX_H