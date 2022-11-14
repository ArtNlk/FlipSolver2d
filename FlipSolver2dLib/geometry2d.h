#ifndef GEOMETRY2D_H
#define GEOMETRY2D_H

#include <vector>
#include <cmath>
#include <algorithm>

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

    friend Vertex operator-(Vertex lhs, Vertex rhs);

    friend Vertex operator*(Vertex lhs, float rhs);

    friend Vertex operator/(Vertex lhs, float rhs);

    friend Vertex operator*(float lhs, Vertex rhs);

    friend Vertex operator+(Vertex lhs, Vertex rhs);

    inline float dot(Vertex rhs)
    {
        return m_x * rhs.m_x + m_y * rhs.m_y + m_z * rhs.m_z;
    }

    inline float distFromZero() {return std::sqrt(m_x * m_x + m_y * m_y + m_z * m_z);}

    //Adapted from https://iquilezles.org/articles/distfunctions2d/
    inline float distToLine(Vertex a, Vertex b)
    {
        Vertex pa = *this-a;
        Vertex len = b-a;
        float h = std::clamp( pa.dot(len)/len.dot(len), 0.0f, 1.0f);
        return (pa - len*h).distFromZero();
    }

    inline Vertex normalized()
    {
        float length = distFromZero();
        if(std::abs(length) < 1e-6f)
        {
            return Vertex(0.f,0.f,0.f);
        }
        return Vertex(m_x / length, m_y / length, m_z / length);
    }

protected:
    float m_x;
    float m_y;
    float m_z;
};

class Geometry2d
{
public:
    Geometry2d();
    Geometry2d(std::vector<Vertex> &verts);

    void addVertex(Vertex v);

    std::vector<Vertex> verts();

    int vertextCount();

    float signedDistance(Vertex point);

    float signedDistance(float x, float y);

protected:
    std::vector<Vertex> m_verts;
};

#endif // GEOMETRY2D_H
