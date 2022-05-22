#include "geometry2d.h"

Geometry2d::Geometry2d()
{
}

Geometry2d::Geometry2d(std::vector<Vertex> &verts)
{
    m_verts = verts;
}

void Geometry2d::addVertex(Vertex v)
{
    m_verts.push_back(v);
}

std::vector<Vertex> Geometry2d::verts()
{
    return m_verts;
}

int Geometry2d::vertextCount()
{
    return m_verts.size();
}

//Adapted from https://iquilezles.org/articles/distfunctions2d/
float Geometry2d::signedDistance(Vertex point)
{
    float sqrdDist = (point-m_verts[0]).dot(point-m_verts[0]);
    float s = 1.0;
    int N = m_verts.size();
    for( int i=0, j=N-1; i<N; j=i, i++ )
    {
        Vertex e = m_verts[j] - m_verts[i];
        Vertex w = point - m_verts[i];
        Vertex b = w - e*std::clamp( w.dot(e)/e.dot(e), 0.0f, 1.0f);
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
    return signedDistance(Vertex(x,y));
}
