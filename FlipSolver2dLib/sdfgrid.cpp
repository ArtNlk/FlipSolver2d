#include "sdfgrid.h"
#include "mathfuncs.h"

SdfGrid::SdfGrid(size_t sizeI, size_t sizeJ) :
    Grid2d(sizeI,sizeJ, 0.f, OOBStrategy::OOB_EXTEND)
{

}

Vec3 SdfGrid::closestSurfacePoint(float i, float j)
{
    return closestSurfacePoint(Vec3(i,j));
}

Vec3 SdfGrid::closestSurfacePoint(Vec3 pos)
{
    Vec3 closestPoint = pos;
    auto* grid = dynamic_cast<Grid2d<float>*>(this);
    float value = simmath::lerpCenteredGrid(pos.x(),pos.y(),*grid);
    Vec3 grad = simmath::gradCenteredGrid(pos.x(),pos.y(),*grid);
    static const int iterationLimit = 100;
    static const int internalIterationLimit = 10;
    for(int i = 0; i < iterationLimit; i++)
    {
        float alpha = 1;
        for(int j = 0; j < internalIterationLimit; j++)
        {
            Vec3 q = closestPoint - alpha*value*grad;
            if(std::abs(simmath::lerpCenteredGrid(q.x(),q.y(),*grid)) < std::abs(value))
            {
                closestPoint = q;
                value = simmath::lerpCenteredGrid(q.x(),q.y(),*grid);
                grad = simmath::gradCenteredGrid(q.x(),q.y(),*grid);
                if(std::abs(value) < 1e-5f)
                {
                    return closestPoint;
                }
            }
            else
            {
                alpha *= 0.7f;
            }
        }
    }
    return closestPoint;
}
