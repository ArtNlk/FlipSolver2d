#include "sdfgrid.h"
#include "mathfuncs.h"

SdfGrid::SdfGrid(int sizeI, int sizeJ) :
    Grid2d(sizeI,sizeJ, 0.f, OOBStrategy::OOB_EXTEND)
{

}

Vertex SdfGrid::closestSurfacePoint(float i, float j)
{
    return closestSurfacePoint(Vertex(i,j));
}

Vertex SdfGrid::closestSurfacePoint(Vertex pos)
{
    Vertex closestPoint = pos;
    auto* grid = dynamic_cast<Grid2d<float>*>(this);
    float value = simmath::lerpCenteredGrid(pos.x(),pos.y(),*grid);
    Vertex grad = simmath::gradCenteredGrid(pos.x(),pos.y(),*grid);
    static const int iterationLimit = 100;
    static const int internalIterationLimit = 10;
    for(int i = 0; i < iterationLimit; i++)
    {
        float alpha = 1;
        for(int j = 0; j < internalIterationLimit; j++)
        {
            Vertex q = closestPoint - alpha*value*grad;
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
