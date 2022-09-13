#include "obstacle.h"

Obstacle::Obstacle(float friction, Geometry2d& geo) :
    m_friction(friction),
    m_geometry(geo)
{
}

float Obstacle::friction()
{
    return m_friction;
}

void Obstacle::setFriction(float newFriction)
{
    m_friction = newFriction;
}

Geometry2d &Obstacle::geometry()
{
    return m_geometry;
}

void Obstacle::setGeometry(Geometry2d &newGeometry)
{
    m_geometry = newGeometry;
}
