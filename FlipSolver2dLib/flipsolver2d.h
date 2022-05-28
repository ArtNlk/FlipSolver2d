#ifndef FLIPSOLVER_H
#define FLIPSOLVER_H

#include <vector>
#include <random>

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

    void updateSdf();

    void updateSolids();

    void updateSources();

    void updateSinks();

    int gridSizeI();

    int gridSizeJ();

    void addGeometry(Geometry2d& geometry);

    void addSource(Geometry2d& geometry);

    void addSink(Geometry2d& geometry);

    void addMarkerParticle(Vertex particle);

    void reseedParticles();

    std::vector<Geometry2d> &geometryObjects();

    std::vector<Geometry2d> &sourceObjects();

    std::vector<Geometry2d> &sinkObjects();

    std::vector<Vertex> &markerParticles();

protected:

    void calcRhs(std::vector<double> &rhs);

    Vertex jitteredPosInCell(int i, int j);

    void countParticles(Grid2d<int> &output);

    int m_extrapolationRadius;
    bool m_useVonNeumannNeighborhood;
    MACFluidGrid m_grid;
    std::vector<Vertex> m_markerParticles;
    PCGSolver m_pcgSolver;
    std::mt19937 m_randEngine;
    std::vector<Geometry2d> m_geometry;
    std::vector<Geometry2d> m_sources;
    std::vector<Geometry2d> m_sinks;
};

#endif // FLIPSOLVER_H
