#ifndef FLIPSOLVER_H
#define FLIPSOLVER_H

#include <vector>

#include "pcgsolver.h"
#include "fluidgrid.h"
#include "uppertriangularmatrix.h"
#include "geometry2d.h"

class FlipSolver
{
public:
    FlipSolver(int sizeX, int sizeY, double fluidDensity, double timestepSize, double sideLength, int extrapRadius = 1, bool vonNeumannNeighbors = false);

    inline MACFluidGrid &grid() {return m_grid;}

    inline void setExrapolationRadius(int radius)
    {
        m_extrapolationRadius = radius;
    }

    void extrapolateVelocityField();

    void project();

    int gridSizeI();

    int gridSizeJ();

    void addGeometry(Geometry2d geometry);

    std::vector<Geometry2d> &geometryObjects();

protected:

    void calcRhs(std::vector<double> &rhs);

    int m_extrapolationRadius;
    bool m_useVonNeumannNeighborhood;
    MACFluidGrid m_grid;
    PCGSolver m_pcgSolver;
    std::vector<Geometry2d> m_geometry;
};

#endif // FLIPSOLVER_H
