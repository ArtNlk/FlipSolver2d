#ifndef FLIPSOLVER_H
#define FLIPSOLVER_H

#include <vector>
#include <random>
#include <memory>

#include "linearindexable2d.h"
#include "materialgrid.h"
#include "obstacle.h"
#include "pcgsolver.h"
#include "fluidgrid.h"
#include "sdfgrid.h"
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

struct PairHash {
public:
  template <typename T, typename U>
  std::size_t operator()(const std::pair<T, U> &x) const
  {
    std::size_t lhs = std::hash<T>()(x.first);
    std::size_t rhs = std::hash<U>()(x.second);
    return rhs + 0x9e3779b9 + (lhs << 6) + (lhs >> 2);
  }
};

class FlipSolver : public LinearIndexable2d
{
public:
    FlipSolver(int sizeI, int sizeJ, int extrapRadius = 1, bool vonNeumannNeighbors = false);

    virtual ~FlipSolver() = default;

    inline void setExrapolationRadius(int radius)
    {
        m_extrapolationRadius = radius;
    }

    int particleCount();

    int cellCount();

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

protected:

    friend class ::LiquidRenderApp;

    virtual double divergenceAt(int i, int j);

    std::vector<int> validSolidNeighborIds(int i, int j);

    void resetGrids();

    virtual void project();

    void applyViscosity();

    virtual void advect();

    virtual void particleUpdate();

    virtual void step();

    virtual void calcPressureRhs(std::vector<double> &rhs);

    void calcViscosityRhs(std::vector<double> &rhs);

    virtual void extrapolateVelocityField(Grid2d<float> &extrapGrid, Grid2d<bool> &flagGrid, int steps = 10);

    virtual DynamicUpperTriangularSparseMatrix getPressureProjectionMatrix();

    DynamicUpperTriangularSparseMatrix getViscosityMatrix();

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

    void updateLinearFluidViscosityMapping();

    void updateValidULinearMapping();

    void updateValidVLinearMapping();

    int linearViscosityVelocitySampleIndexU(int i, int j);

    int linearViscosityVelocitySampleIndexV(int i, int j);

    float maxParticleVelocity();

    int m_extrapolationRadius;
    bool m_useVonNeumannNeighborhood;
    int m_frameNumber;
    int m_validVVelocitySampleCount;
    int m_validUVelocitySampleCount;
    std::vector<MarkerParticle> m_markerParticles;
    PCGSolver m_pcgSolver;
    std::mt19937 m_randEngine;
    std::vector<Obstacle> m_obstacles;
    std::vector<Emitter> m_sources;
    std::vector<Sink> m_sinks;
    std::vector<Emitter> m_initialFluid;

    StaggeredVelocityGrid m_fluidVelocityGrid;
    StaggeredVelocityGrid m_savedFluidVelocityGrid;
    MaterialGrid m_materialGrid;
    SdfGrid m_solidSdf;
    SdfGrid m_fluidSdf;
    Grid2d<bool> m_knownCenteredParams;
    Grid2d<float> m_viscosityGrid;
    Grid2d<int> m_emitterId;
    Grid2d<int> m_solidId;
    Grid2d<int> m_fluidParticleCounts;
    Grid2d<float> m_divergenceControl;
    Grid2d<float> m_testGrid;

    std::unordered_map<std::pair<int,int>,int,PairHash> m_uVelocitySamplesMap;
    std::unordered_map<std::pair<int,int>,int,PairHash> m_vVelocitySamplesMap;
};

#endif // FLIPSOLVER_H
