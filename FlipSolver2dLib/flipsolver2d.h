#ifndef FLIPSOLVER_H
#define FLIPSOLVER_H

#include <vector>
#include <random>
#include <memory>

#include "obstacle.h"
#include "pcgsolver.h"
#include "fluidgrid.h"
#include "uppertriangularmatrix.h"
#include "geometry2d.h"
#include "emitter.h"

class LiquidRenderApp;

struct MarkerParticle
{
    Vertex position;
    Vertex velocity;
    float viscosity;
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

    void init();

    void extrapolateVelocityField(int steps = std::numeric_limits<int>().max());

    UpperTriangularMatrix getPressureProjectionMatrix();

    UpperTriangularMatrix getViscosityMatrix();

    void project();

    void applyViscosity();

    void advect();

    void particleUpdate(Grid2d<float>& prevU, Grid2d<float>& prevV);

    void step();

    void stagedStep();

    void stepFrame();

    void updateSolids();

    void updateSources();

    void updateSinks();

    void updateInitialFluid();

    int gridSizeI();

    int gridSizeJ();

    void addGeometry(Obstacle &geometry);

    void addSource(Emitter& emitter);

    void addSink(Geometry2d& geometry);

    void addInitialFluid(Emitter& emitter);

    void addMarkerParticle(Vertex particle);

    void addMarkerParticle(MarkerParticle particle);

    int frameNumber();

    std::vector<Obstacle> &geometryObjects();

    std::vector<Emitter> &sourceObjects();

    std::vector<Geometry2d> &sinkObjects();

    std::vector<MarkerParticle> &markerParticles();

    SimulationStepStage stepStage();

protected:

    friend class ::LiquidRenderApp;

    void calcPressureRhs(std::vector<double> &rhs);

    void calcViscosityRhs(std::vector<double> &rhs);

    Vertex jitteredPosInCell(int i, int j);

    void reseedParticles(Grid2d<int> &particleCounts);

    void seedInitialFluid();

    void countParticles(Grid2d<int> &output);

    void updateMaterialsFromParticles(Grid2d<int> &particleCount);

    void updateVelocityFromSolids();

    void applyPressuresToVelocityField(std::vector<double> pressures);

    Vertex rk3Integrate(Vertex currentPosition, float dt);

    void particleToGrid();

    void applyGlobalAcceleration();

    float maxParticleVelocity();

    int m_extrapolationRadius;
    bool m_useVonNeumannNeighborhood;
    int m_frameNumber;
    MACFluidGrid m_grid;
    std::vector<MarkerParticle> m_markerParticles;
    PCGSolver m_pcgSolver;
    std::mt19937 m_randEngine;
    std::vector<Obstacle> m_obstacles;
    std::vector<Emitter> m_sources;
    std::vector<Geometry2d> m_sinks;
    std::vector<Emitter> m_initialFluid;
    SimulationStepStage m_stepStage;
};

#endif // FLIPSOLVER_H
