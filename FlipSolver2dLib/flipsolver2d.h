#ifndef FLIPSOLVER_H
#define FLIPSOLVER_H

#include <array>
#include <chrono>
#include <vector>
#include <random>
#include <memory>

#include "PressureIPPCoeficients.h"
#include "emitter.h"
#include "geometry2d.h"
#include "inversepoissonpreconditioner.h"
#include "linearindexable2d.h"
#include "linearsolver.h"
#include "markerparticlesystem.h"
#include "materialgrid.h"
#include "obstacle.h"
#include "pressuredata.h"
#include "sdfgrid.h"
#include "sink.h"
#include "staggeredvelocitygrid.h"
#include "staticmatrix.h"
#include "threadpool.h"
#include "viscositymodel.h"

#include <Eigen/Sparse>

enum SimulationMethod : char {SIMULATION_LIQUID, SIMULATION_SMOKE, SIMULATION_FIRE, SIMULATION_NBFLIP};

enum ParameterHandlingMethod : char {PARTICLE, HYBRID, GRID};

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
    bool useHeavyViscosity;
    SimulationMethod simulationMethod;
    ParameterHandlingMethod parameterHandlingMethod;
};

enum SolverStage {
    ADVECTION = 0,
    DECOMPOSITION,
    DENSITY,
    PARTICLE_REBIN,
    PARTICLE_TO_GRID,
    GRID_UPDATE,
    AFTER_TRANSFER,
    PRESSURE,
    VISCOSITY,
    REPRESSURE,
    PARTICLE_UPDATE,
    PARTICLE_RESEED,
    SOLVER_STAGE_COUNT
};

template<class IterableEnum, IterableEnum endVal>
inline IterableEnum& incrementEnum(IterableEnum& state) {
    const int i = static_cast<int>(state)+1;
    state = static_cast<IterableEnum>((i) % endVal);
    return state;
}

template<class IterableEnum, IterableEnum endVal>
inline bool nextEnum(IterableEnum& state) {
    bool isLast = state == (endVal-1);
    state = incrementEnum<IterableEnum, endVal>(state);
    return !isLast;
}

class SolverTimeStats {
public:
    using Clock = std::chrono::high_resolution_clock;
    using TimePoint = std::chrono::time_point<Clock,Clock::duration>;
    using MsDuration = std::chrono::duration<float, std::milli>;
    using StageTimings = std::array<float, SOLVER_STAGE_COUNT>;

    SolverTimeStats()
    {
        reset();
    }

    void reset(){
        m_times.fill(0.f);
        m_lastTimePoint = Clock::now();
        m_frameStartTimePoint = Clock::now();
        m_substepsTaken = 0;
    }

    void endStage(SolverStage s){
        TimePoint currentTime = Clock::now();
        m_times.at(s) += MsDuration(currentTime - m_lastTimePoint).count();
        m_lastTimePoint = Clock::now();
    }

    void endFrame()
    {
        TimePoint currentTime = Clock::now();
        m_totalFrameTime = MsDuration(currentTime - m_frameStartTimePoint).count();
    }

    void addSubstep()
    {
        m_substepsTaken++;
    }

    int substepCount() const
    {
        return m_substepsTaken;
    }

    StageTimings timings() const
    {
        return m_times;
    }

    float frameTime() const
    {
        return m_totalFrameTime;
    }

protected:
    StageTimings m_times;
    TimePoint m_lastTimePoint;
    TimePoint m_frameStartTimePoint;
    int m_substepsTaken;
    float m_totalFrameTime;
};

class FlipSolver : public LinearIndexable2d
{
public:
    FlipSolver(const FlipSolverParameters *p);

    virtual ~FlipSolver();

    size_t particleCount();

    size_t cellCount();

    void stepFrame();

    void updateSolids();

    void updateSources();

    void updateSinks();

    void updateInitialFluid();

    size_t gridSizeI();

    size_t gridSizeJ();

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

    const SolverTimeStats &timeStats() const;

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

    void eulerAdvectionThread(Range range,
                              const Grid2d<float> &inputGrid,
                              Grid2d<float> &outputGrid);

protected:
    void pruneParticles();

    virtual double divergenceAt(size_t i, size_t j);

    std::vector<int> validSolidNeighborIds(ssize_t i, ssize_t j);

    virtual void firstFrameInit();

    virtual void project();

    void applyViscosity();

    virtual void advect();

    virtual void eulerAdvectParameters(){};

    void densityCorrection();

    void updateDensityGrid();

    void updateDensityGridThread(Range r, Grid2d<float>& centeredWeights);

    void adjustParticlesByDensity();

    void adjustParticlesByDensityThread(Range r);

    void advectThread(Range range, std::vector<std::vector<RebinRecord>> &rebinningSets);

    virtual void particleUpdate();

    virtual void afterTransfer();

    virtual void step();

    virtual void calcPressureRhs(std::vector<double> &rhs);

    void calcDensityCorrectionRhs(std::vector<double> &rhs);
    
    virtual IndexedPressureParameters getPressureProjectionMatrix();

    virtual IndexedIPPCoefficients getIPPCoefficients(const IndexedPressureParameters& mat);

    Vertex jitteredPosInCell(size_t i, size_t j);

    virtual void reseedParticles();

    virtual void seedInitialFluid();

    virtual void countParticles();

    virtual void updateMaterials();

    void updateVelocityFromSolids();

    virtual void applyPressuresToVelocityField(const std::vector<double> &pressures);

    void applyPressureThreadU(Range range, const std::vector<double> &pressures);

    void applyPressureThreadV(Range range,const std::vector<double> &pressures);
    
    Vertex rk4Integrate(Vertex currentPosition, StaggeredVelocityGrid &grid, float dt);

    virtual void gridUpdate();

    virtual void particleToGrid();

    virtual void applyBodyForces();

    virtual void updateSdf();

    void updateSdfThread(Range range);

    virtual void particleVelocityToGrid();

    void particleVelocityToGridThread(Range r);

    virtual void centeredParamsToGrid();

    void centeredParamsToGridThread(Range r, Grid2d<float>& cWeights);

    void extrapolateLevelsetInside(SdfGrid& grid);

    void extrapolateLevelsetOutside(SdfGrid& grid);

    float maxParticleVelocity();

    float maxGridVelocity();

    int m_frameNumber;
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
    Eigen::VectorXd m_rhs;
    Eigen::VectorXd m_solverResult;
    std::vector<double> m_pressureRhs;
    std::vector<double> m_pressureSolverResult;

    std::shared_ptr<ViscosityModel> m_viscosityModel;

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
    ParameterHandlingMethod m_parameterHandlingMethod;
    SolverTimeStats m_stats;

    //using precond = Eigen::IncompleteCholesky<double,Eigen::Upper>;
    using precond = InversePoissonPreconditioner<double, Eigen::Upper>;
    IndexedPressureParameters m_pressureMatrix;
    IndexedIPPCoefficients m_pressurePrecond;
    // Eigen::ConjugateGradient<Eigen::SparseMatrix<double>,Eigen::Upper,precond> m_pressureSolver;
    LinearSolver m_pressureSolver;

    size_t m_testValuePropertyIndex;
    size_t m_viscosityPropertyIndex;

    float m_frameTime;
    float m_avgFrameMs;
};

#endif // FLIPSOLVER_H
