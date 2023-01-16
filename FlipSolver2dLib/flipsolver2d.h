#ifndef FLIPSOLVER_H
#define FLIPSOLVER_H

#include <vector>
#include <random>
#include <memory>

#include "fluidcell.h"
#include "obstacle.h"
#include "pcgsolver.h"
#include "fluidgrid.h"
#include "uppertriangularmatrix.h"
#include "geometry2d.h"
#include "emitter.h"
#include "sink.h"

class LiquidRenderApp;

struct MarkerParticle
{
    Vertex position;
    Vertex velocity;
    float viscosity;
    float temperature;
    float smokeConcentrartion;
    FluidMaterial material = FluidMaterial::FLUID;
    float testValue = 0.f;
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

    virtual ~FlipSolver() = default;

    inline MACFluidGrid &grid() {return m_grid;}

    inline void setExrapolationRadius(int radius)
    {
        m_extrapolationRadius = radius;
    }

    void init();

    virtual void extrapolateVelocityField(Grid2d<float> &extrapGrid, Grid2d<bool> &flagGrid, int steps = 10);

    virtual DynamicUpperTriangularSparseMatrix getPressureProjectionMatrix();

    DynamicUpperTriangularSparseMatrix getViscosityMatrix();

    int particleCount();

    virtual void project();

    void applyViscosity();

    virtual void advect();

    virtual void particleUpdate();

    virtual void step();

    void stepFrame();

    void updateSolids();

    void updateSources();

    void updateSinks();

    void updateInitialFluid();

    int gridSizeI();

    int gridSizeJ();

    void addGeometry(Obstacle &geometry);

    void addSource(Emitter& emitter);

    void addSink(Sink &geometry);

    void addInitialFluid(Emitter& emitter);

    void addMarkerParticle(Vertex particle);

    void addMarkerParticle(MarkerParticle particle);

    int frameNumber();

    std::vector<Obstacle> &geometryObjects();

    std::vector<Emitter> &sourceObjects();

    std::vector<Sink> &sinkObjects();

    std::vector<MarkerParticle> &markerParticles();

    SimulationStepStage stepStage();

protected:

    friend class ::LiquidRenderApp;

    virtual void calcPressureRhs(std::vector<double> &rhs);

    void calcViscosityRhs(std::vector<double> &rhs);

    Vertex jitteredPosInCell(int i, int j);

    virtual void reseedParticles();

    virtual void seedInitialFluid();

    virtual void countParticles();

    virtual void updateMaterials();

    void updateVelocityFromSolids();

    virtual void applyPressuresToVelocityField(std::vector<double> &pressures);

    Vertex rk3Integrate(Vertex currentPosition, float dt, StaggeredVelocityGrid &grid);

    virtual void particleToGrid();

    virtual void applyBodyForces();

    virtual void updateSdf();

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
    std::vector<Sink> m_sinks;
    std::vector<Emitter> m_initialFluid;
    SimulationStepStage m_stepStage;
};

#endif // FLIPSOLVER_H
