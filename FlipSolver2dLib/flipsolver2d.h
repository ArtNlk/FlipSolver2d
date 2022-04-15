#ifndef FLIPSOLVER_H
#define FLIPSOLVER_H

#include "pcgsolver.h"
#include "fluidgrid.h"

class FlipSolver
{
public:
    FlipSolver(int sizeX, int sizeY, double fluidDensity, double timestepSize, double sideLength);

    inline MACFluidGrid &grid() {return m_grid;}

protected:
    void project();

    void calcRhs(std::vector<double> &rhs);

    MACFluidGrid m_grid;
    PCGSolver m_pcgSolver;
};

#endif // FLIPSOLVER_H
