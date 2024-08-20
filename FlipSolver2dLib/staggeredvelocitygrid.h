#ifndef STAGGEREDVELOCITYGRID_H
#define STAGGEREDVELOCITYGRID_H

#include "linearindexable2d.h"
#include "grid2d.h"
#include "geometry2d.h"

class StaggeredVelocityGrid : public LinearIndexable2d
{
public:
    StaggeredVelocityGrid(size_t sizeI, size_t sizeJ);

    Grid2d<float>& velocityGridU();
    Grid2d<float>& velocityGridV();

    const Grid2d<float>& velocityGridU() const;
    const Grid2d<float>& velocityGridV() const;

    Grid2d<bool>& uSampleValidityGrid();
    Grid2d<bool>& vSampleValidityGrid();

    float &u(ssize_t i,ssize_t j);

    float &v(ssize_t i,ssize_t j);

    float getU(ssize_t i,ssize_t j) const;

    float getV(ssize_t i,ssize_t j) const;

    void setU(ssize_t i,ssize_t j, float u);

    void setV(ssize_t i,ssize_t j, float v);

    bool getUValidity(ssize_t i,ssize_t j);

    bool getVValidity(ssize_t i,ssize_t j);

    void setUValidity(ssize_t i,ssize_t j, bool uValidity);

    void setVValidity(ssize_t i,ssize_t j, bool vValidity);

    float &u(Index2d idx);

    float &v(Index2d idx);

    float getU(Index2d idx) const;

    float getV(Index2d idx) const;

    void setU(Index2d idx, float u);

    void setV(Index2d idx, float v);

    bool getUValidity(Index2d idx);

    bool getVValidity(Index2d idx);

    void setUValidity(Index2d idx, bool uValidity);

    void setVValidity(Index2d idx, bool vValidity);

    void extrapolate(int extrapolationRadius);

    Vec3 velocityAt(float i, float j) const;
    Vec3 velocityAt(Vec3 position) const;

protected:
    Grid2d<float> m_velocityGridU;
    Grid2d<float> m_velocityGridV;
    Grid2d<bool> m_uSampleValidity;
    Grid2d<bool> m_vSampleValidity;

    int m_extrapolationNeighborRadius = 1;
};

#endif // STAGGEREDVELOCITYGRID_H
