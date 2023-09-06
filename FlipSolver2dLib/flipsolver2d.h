#ifndef FLIPSOLVER_H
#define FLIPSOLVER_H

#include <vector>
#include <random>
#include <memory>

#include "linearindexable2d.h"
#include "markerparticlesystem.h"
#include "materialgrid.h"
#include "obstacle.h"
#include "linearsolver.h"
#include "sdfgrid.h"
#include "staggeredvelocitygrid.h"
#include "threadpool.h"
#include "uppertriangularmatrix.h"
#include "geometry2d.h"
#include "emitter.h"
#include "sink.h"

enum SimulationMethod : char {SIMULATION_LIQUID, SIMULATION_SMOKE, SIMULATION_FIRE, SIMULATION_NBFLIP};

struct FlipSolverParameters
{
    double fluidDensity;
    unsigned int seed;
    double dx;
    int particlesPerCell;
    Vertex globalAcceleration;
    float resolution;
    int fps;
    int maxSubsteps;
    float picRatio;
    float cflNumber;
    float particleScale;
    int pcgIterLimit;
    float domainSizeI;
    float domainSizeJ;
    int gridSizeI;
    int gridSizeJ;
    float sceneScale;
    bool viscosityEnabled;
    SimulationMethod simulationMethod;
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
    FlipSolver(const FlipSolverParameters *p);

    virtual ~FlipSolver();

    size_t particleCount();

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

    int frameNumber();

    std::vector<Obstacle> &geometryObjects();

    std::vector<Emitter> &sourceObjects();

    std::vector<Sink> &sinkObjects();

    MarkerParticleSystem &markerParticles();

    const MaterialGrid &materialGrid() const;

    const StaggeredVelocityGrid &fluidVelocityGrid() const;

    const SdfGrid &fluidSdf() const;

    const SdfGrid &solidSdf() const;

    const Grid2d<float> &testGrid() const;

    float stepDt() const;

    double dx() const;

    SimulationMethod simulationMethod() const;

    float domainSizeI() const;

    float domainSizeJ() const;

    double fluidDensity() const;

    int particlesPerCell() const;

    float picRatio() const;

    float cflNumber() const;

    int maxSubsteps() const;

    int fps() const;

    float frameDt() const;

    Vertex globalAcceleration() const;

    float sceneScale() const;

    float lastFrameTime() const;

    float avgFrameTime() const;

    size_t testValuePropertyIndex();

    virtual void initAdditionalParameters();

protected:

    virtual double divergenceAt(int i, int j);

    std::vector<int> validSolidNeighborIds(int i, int j);

    virtual void firstFrameInit();

    virtual void project();

    virtual LinearSolver::MatElementProvider getPressureMatrixElementProvider();

    LinearSolver::SparseMatRowElements getMatFreeElementForLinIdx(unsigned int i);

    void applyViscosity();

    virtual void advect();

    void densityCorrection();

    void updateDensityGrid();

    void updateDensityGridThread(Range r, Grid2d<float>& centeredWeights);

    void adjustParticlesByDensity();

    void adjustParticlesByDensityThread(Range r);

    void advectThread(Range range);

    virtual void particleUpdate();

    virtual void afterTransfer();

    virtual void step();

    virtual void calcPressureRhs(std::vector<double> &rhs);

    void calcViscosityRhs(std::vector<double> &rhs);

    void calcDensityCorrectionRhs(std::vector<double> &rhs);

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

    void applyPressureThreadU(Range range, const std::vector<double> &pressures);

    void applyPressureThreadV(Range range,const std::vector<double> &pressures);
    
    Vertex rk4Integrate(Vertex currentPosition, StaggeredVelocityGrid &grid);

    virtual void particleToGrid();

    virtual void applyBodyForces();

    virtual void updateSdf();

    void updateSdfThread(Range range);

    virtual void particleVelocityToGrid();

    void particleVelocityToGridThread(Range r, Grid2d<float>& uWeights, Grid2d<float>& vWeights);

    virtual void centeredParamsToGrid();

    void extrapolateLevelsetInside(SdfGrid& grid);

    void extrapolateLevelsetOutside(SdfGrid& grid);

    float maxParticleVelocity();

    float maxGridVelocity();

    int m_frameNumber;
    int m_validVVelocitySampleCount;
    int m_validUVelocitySampleCount;
    MarkerParticleSystem m_markerParticles;
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
    Grid2d<float> m_densityGrid;
    Grid2d<float> m_testGrid;
    std::vector<double> m_rhs;
    Grid2d<double> m_pressures;

    LinearSolver m_pcgSolver;
    std::shared_ptr<IPreconditioner> m_projectPreconditioner;
    std::shared_ptr<IPreconditioner> m_densityPreconditioner;
    std::shared_ptr<IPreconditioner> m_viscosityPreconditioner;

    float m_stepDt;
    float m_frameDt;
    double m_dx;
    double m_fluidDensity;
    unsigned int m_seed;
    int m_particlesPerCell;
    Vertex m_globalAcceleration;
    float m_resolution;
    int m_fps;
    int m_maxSubsteps;
    float m_picRatio;
    float m_cflNumber;
    float m_particleScale;
    int m_pcgIterLimit;
    float m_domainSizeI;
    float m_domainSizeJ;
    float m_sceneScale;
    double m_projectTolerance;
    bool m_viscosityEnabled;
    SimulationMethod m_simulationMethod;

    size_t m_testValuePropertyIndex;
    size_t m_viscosityPropertyIndex;

    float m_frameTime;
    float m_avgFrameMs;

    std::unordered_map<std::pair<int,int>,int,PairHash> m_uVelocitySamplesMap;
    std::unordered_map<std::pair<int,int>,int,PairHash> m_vVelocitySamplesMap;
};

#endif // FLIPSOLVER_H
