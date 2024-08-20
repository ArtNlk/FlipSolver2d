#include "geometry2d.h"

Vec3 operator-(Vec3 lhs, Vec3 rhs)
{
    return Vec3(lhs.m_x - rhs.m_x, lhs.m_y - rhs.m_y, lhs.m_z - rhs.m_z);
}

Vec3 operator*(Vec3 lhs, float rhs)
{
    return Vec3(lhs.m_x*rhs, lhs.m_y*rhs, lhs.m_z*rhs);
}

Vec3 operator/(Vec3 lhs, float rhs)
{
    return Vec3(lhs.m_x/rhs, lhs.m_y/rhs, lhs.m_z/rhs);
}


Vec3 operator*(float lhs, Vec3 rhs)
{
    return Vec3(rhs.m_x*lhs, rhs.m_y*lhs, rhs.m_z*lhs);
}

Vec3 operator+(Vec3 lhs, Vec3 rhs)
{
    return Vec3(lhs.m_x + rhs.m_x, lhs.m_y + rhs.m_y, lhs.m_z + rhs.m_z);
}

Geometry2d::Geometry2d()
{
}

Geometry2d::Geometry2d(std::vector<Vec3> &verts)
{
    m_verts = verts;
}

void Geometry2d::addVertex(Vec3 v)
{
    m_verts.push_back(v);
}

std::vector<Vec3> Geometry2d::verts()
{
    return m_verts;
}

int Geometry2d::vertextCount()
{
    return m_verts.size();
}

//Adapted from https://iquilezles.org/articles/distfunctions2d/
float Geometry2d::signedDistance(Vec3 point)
{
    float sqrdDist = (point-m_verts[0]).dot(point-m_verts[0]);
    float s = 1.0;
    size_t N = m_verts.size();
    for(size_t i=0, j=N-1; i<N; j=i, i++ )
    {
        Vec3 e = m_verts[j] - m_verts[i];
        Vec3 w = point - m_verts[i];
        Vec3 b = w - e*std::clamp( w.dot(e)/e.dot(e), 0.0f, 1.0f);
        sqrdDist = std::min( sqrdDist, b.dot(b) );
        bool xf = point.y()>=m_verts[i].y();
        bool yf = point.y()<m_verts[j].y();
        bool zf = e.x()*w.y()>e.y()*w.x();
        bool allC = xf && yf && zf;
        bool allNC = !xf && !yf && !zf;
        if( allC || allNC ) s*=-1.0;
    }
    return s*sqrt(sqrdDist);
}

float Geometry2d::signedDistance(float x, float y)
{
    return signedDistance(Vec3(x,y));
}

Vec3::Vec3(float x, float y, float z) :
    m_x(x),
    m_y(y),
    m_z(z)
{
}

Vec3::Vec3(std::pair<float, float> &p) :
    m_x(p.first),
    m_y(p.second),
    m_z(0.f)
{

}
