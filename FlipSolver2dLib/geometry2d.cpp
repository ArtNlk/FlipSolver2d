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
