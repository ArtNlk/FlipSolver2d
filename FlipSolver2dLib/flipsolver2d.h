#ifndef FLIPSOLVER_H
#define FLIPSOLVER_H

#include <vector>
#include <random>

#include "pcgsolver.h"
#include "fluidgrid.h"
#include "uppertriangularmatrix.h"
#include "geometry2d.h"

class LiquidRenderApp;

struct MarkerParticle
{
    Vertex position;
    Vertex velocity;
};

enum SimulationStepStage : int {STAGE_RESEED,
                                STAGE_UPDATE_MATERIALS,
                                STAGE_TRANSFER_PARTICLE_VELOCITY,
                                STAGE_EXTRAPOLATE_VELOCITY,
                                STAGE_APPLY_GLOBAL_ACCELERATION,
                                STAGE_PROJECT,
                                STAGE_EXTRAPOLATE_AFTER_PROJECTION,
                                STAGE_UPDATE_PARTICLE_VELOCITIES,
                                STAGE_ADVECT,
                                STAGE_ITER_END};
inline SimulationStepStage& operator++(SimulationStepStage& state, int) {
    const int i = static_cast<int>(state)+1;
    state = static_cast<SimulationStepStage>((i) % STAGE_ITER_END);
    return state;
}

class FlipSolver
{
public:
    FlipSolver(int extrapRadius = 1, bool vonNeumannNeighbors = false);

    inline MACFluidGrid &grid() {return m_grid;}

    inline void setExrapolationRadius(int radius)
    {
        m_extrapolationRadius = radius;
    }

    void extrapolateVelocityField(int steps = std::numeric_limits<int>().max());

    void project();

    void advect();

    void step();

    void stagedStep();

    void stepFrame();

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

    void addMarkerParticle(MarkerParticle particle);

    std::vector<Geometry2d> &geometryObjects();

    std::vector<Geometry2d> &sourceObjects();

    std::vector<Geometry2d> &sinkObjects();

    std::vector<MarkerParticle> &markerParticles();

    SimulationStepStage stepStage();

protected:

    friend class ::LiquidRenderApp;

    void calcRhs(std::vector<double> &rhs);

    Vertex jitteredPosInCell(int i, int j);

    void reseedParticles(Grid2d<int> &particleCounts);

    void countParticles(Grid2d<int> &output);

    void updateMaterialsFromParticles(Grid2d<int> &particleCount);

    Vertex rk3Integrate(Vertex currentPosition, float dt);

    void particleVelocityToGrid(Grid2d<float> &gridU, Grid2d<float> &gridV);

    void applyGlobalAcceleration();

    int m_extrapolationRadius;
    bool m_useVonNeumannNeighborhood;
    MACFluidGrid m_grid;
    std::vector<MarkerParticle> m_markerParticles;
    PCGSolver m_pcgSolver;
    std::mt19937 m_randEngine;
    std::vector<Geometry2d> m_geometry;
    std::vector<Geometry2d> m_sources;
    std::vector<Geometry2d> m_sinks;
    SimulationStepStage m_stepStage;
};

#endif // FLIPSOLVER_H
