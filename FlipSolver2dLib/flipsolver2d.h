#ifndef FLIPSOLVER_H
#define FLIPSOLVER_H

#include "pcgsolver.h"
#include "fluidgrid.h"
#include "uppertriangularmatrix.h"

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

protected:

    void calcRhs(std::vector<double> &rhs);

    int m_extrapolationRadius;
    bool m_useVonNeumannNeighborhood;
    MACFluidGrid m_grid;
    PCGSolver m_pcgSolver;
};

#endif // FLIPSOLVER_H
