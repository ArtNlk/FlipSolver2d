#ifndef FLIPSOLVER_H
#define FLIPSOLVER_H

#include <array>
#include <chrono>
#include <vector>
#include <random>
#include <memory>

#include "inversepoissonpreconditioner.h"
#include "linearindexable2d.h"
#include "markerparticlesystem.h"
#include "materialgrid.h"
#include "obstacle.h"
#include "linearsolver.h"
#include "sdfgrid.h"
#include "staggeredvelocitygrid.h"
#include "threadpool.h"
#include "staticmatrix.h"
#include "geometry2d.h"
#include "emitter.h"
#include "sink.h"

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
    SimulationMethod simulationMethod;
    ParameterHandlingMethod parameterHandlingMethod;
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

    virtual double divergenceAt(int i, int j);

    std::vector<int> validSolidNeighborIds(int i, int j);

    virtual void firstFrameInit();

    virtual void project();

    virtual LinearSolver::MatElementProvider getPressureMatrixElementProvider();

    LinearSolver::SparseMatRowElements getMatFreeElementForLinIdx(unsigned int i);

    void applyViscosity();

    virtual void advect();

    virtual void eulerAdvectParameters(){};

    void densityCorrection();

    void updateDensityGrid();

    void updateDensityGridThread(Range r, Grid2d<float>& centeredWeights);

    void adjustParticlesByDensity();

    void adjustParticlesByDensityThread(Range r);

    void advectThread(Range range, std::vector<std::vector<RebinRecord> > &rebinningSets);

    virtual void particleUpdate();

    virtual void afterTransfer();

    virtual void step();

    virtual void calcPressureRhs(Eigen::VectorXd &rhs);

    void calcViscosityRhs(Eigen::VectorXd &rhs, Grid2d<float> &sourceGrid);

    void calcDensityCorrectionRhs(Eigen::VectorXd &rhs);
    
    virtual Eigen::SparseMatrix<double, Eigen::RowMajor> getPressureProjectionMatrix();
    
    Eigen::SparseMatrix<double, Eigen::RowMajor> getViscosityMatrix();

    Vertex jitteredPosInCell(int i, int j);

    virtual void reseedParticles();

    virtual void seedInitialFluid();

    virtual void countParticles();

    virtual void updateMaterials();

    void updateVelocityFromSolids();

    virtual void applyPressuresToVelocityField(Eigen::VectorXd &pressures);

    void applyPressureThreadU(Range range, const Eigen::VectorXd &pressures);

    void applyPressureThreadV(Range range,const Eigen::VectorXd &pressures);
    
    Vertex rk4Integrate(Vertex currentPosition, StaggeredVelocityGrid &grid, float dt);

    virtual void gridUpdate();

    virtual void particleToGrid();

    virtual void applyBodyForces();

    virtual void updateSdf();

    void updateSdfThread(Range range);

    virtual void particleVelocityToGrid();

    void particleVelocityToGridThread(Range r, Grid2d<float>& uWeights, Grid2d<float>& vWeights);

    virtual void centeredParamsToGrid();

    void centeredParamsToGridThread(Range r, Grid2d<float>& cWeights);

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
    Eigen::VectorXd m_rhs;
    Eigen::VectorXd m_solverResult;

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
    ParameterHandlingMethod m_parameterHandlingMethod;
    SolverTimeStats m_stats;

    //using precond = Eigen::IncompleteCholesky<double,Eigen::Upper>;
    using precond = InversePoissonPreconditioner<double, Eigen::Upper>;
    Eigen::SparseMatrix<double,Eigen::RowMajor> m_pressureMatrix;
    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>,Eigen::Upper,precond> m_pressureSolver;

    //Eigen::SparseMatrix<double,Eigen::RowMajor> m_viscosityMatrix;
    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Upper> m_viscositySolver;

    size_t m_testValuePropertyIndex;
    size_t m_viscosityPropertyIndex;

    float m_frameTime;
    float m_avgFrameMs;

    std::unordered_map<std::pair<int,int>,int,PairHash> m_uVelocitySamplesMap;
    std::unordered_map<std::pair<int,int>,int,PairHash> m_vVelocitySamplesMap;
};

#endif // FLIPSOLVER_H
