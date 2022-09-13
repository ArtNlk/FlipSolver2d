#ifndef OBSTACLE_H
#define OBSTACLE_H

#include "geometry2d.h"

class Obstacle
{
public:
    Obstacle(float friction, Geometry2d& geo);

    float friction();
    void setFriction(float newFriction);

    Geometry2d &geometry();
    void setGeometry(Geometry2d &newGeometry);

protected:
    float m_friction;
    Geometry2d m_geometry;
};

#endif // OBSTACLE_H
