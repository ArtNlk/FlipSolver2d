#ifndef STAGGEREDVELOCITYGRID_H
#define STAGGEREDVELOCITYGRID_H

#include "linearindexable2d.h"
#include "grid2d.h"
#include "geometry2d.h"

class StaggeredVelocityGrid : public LinearIndexable2d
{
public:
    StaggeredVelocityGrid(int sizeI, int sizeJ);

    Grid2d<float>& velocityGridU();
    Grid2d<float>& velocityGridV();

    Grid2d<bool>& uSampleValidityGrid();
    Grid2d<bool>& vSampleValidityGrid();

    float &u(int i,int j);

    float &v(int i,int j);

    float getU(int i,int j) const;

    float getV(int i,int j) const;

    void setU(int i,int j, float u);

    void setV(int i,int j, float v);

    bool getUValidity(int i,int j);

    bool getVValidity(int i,int j);

    void setUValidity(int i,int j, bool uValidity);

    void setVValidity(int i,int j, bool vValidity);

    void extrapolate(int extrapolationRadius);

    Vertex velocityAt(float i, float j) const;
    Vertex velocityAt(Vertex position) const;

protected:
    Grid2d<float> m_velocityGridU;
    Grid2d<float> m_velocityGridV;
    Grid2d<bool> m_uSampleValidity;
    Grid2d<bool> m_vSampleValidity;

    int m_extrapolationNeighborRadius = 1;
};

#endif // STAGGEREDVELOCITYGRID_H
