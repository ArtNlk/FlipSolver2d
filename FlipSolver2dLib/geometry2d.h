#ifndef GEOMETRY2D_H
#define GEOMETRY2D_H

#include <vector>
#include <cmath>
#include <algorithm>

class Vec3
{
public:
    Vec3(float x = 0.0f, float y = 0.0f, float z = 0.0f);

    Vec3(std::pair<float,float> &p);

    inline float &x() {return m_x;}
    inline float &y() {return m_y;}
    inline float &z() {return m_z;}

    inline const float &x() const {return m_x;}
    inline const float &y() const {return m_y;}
    inline const float &z() const {return m_z;}

    friend Vec3 operator-(Vec3 lhs, Vec3 rhs);

    friend Vec3 operator*(Vec3 lhs, float rhs);

    friend Vec3 operator/(Vec3 lhs, float rhs);

    friend Vec3 operator*(float lhs, Vec3 rhs);

    friend Vec3 operator+(Vec3 lhs, Vec3 rhs);

    inline float dot(Vec3 rhs)
    {
        return m_x * rhs.m_x + m_y * rhs.m_y + m_z * rhs.m_z;
    }

    inline float distFromZero() {return std::sqrt(m_x * m_x + m_y * m_y + m_z * m_z);}

    //Adapted from https://iquilezles.org/articles/distfunctions2d/
    inline float distToLine(Vec3 a, Vec3 b)
    {
        Vec3 pa = *this-a;
        Vec3 len = b-a;
        float h = std::clamp( pa.dot(len)/len.dot(len), 0.0f, 1.0f);
        return (pa - len*h).distFromZero();
    }

    inline Vec3 normalized()
    {
        float length = distFromZero();
        if(std::abs(length) < 1e-6f)
        {
            return Vec3(0.f,0.f,0.f);
        }
        return Vec3(m_x / length, m_y / length, m_z / length);
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
    Geometry2d(std::vector<Vec3> &verts);

    void addVertex(Vec3 v);

    std::vector<Vec3> verts();

    int vertextCount();

    float signedDistance(Vec3 point);

    float signedDistance(float x, float y);

protected:
    std::vector<Vec3> m_verts;
};

#endif // GEOMETRY2D_H
