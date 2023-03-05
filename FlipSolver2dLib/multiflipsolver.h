#ifndef MULTIFLIPSOLVER_H
#define MULTIFLIPSOLVER_H

#include "flipsolver2d.h"
#include "staggeredvelocitygrid.h"

class MultiflipSolver : public FlipSolver
{
public:
    MultiflipSolver(int sizeI, int sizeJ, int extrapRadius = 1, bool vonNeumannNeighbors = false);

    void calcPressureRhs(std::vector<double> &rhs) override;

    DynamicUpperTriangularSparseMatrix getPressureProjectionMatrix() override;

    void step() override;

    void advect() override;

    void updateSdf() override;

    void combineSdf();

    void extrapolateSdf(Grid2d<float>& sdfGrid);

    void particleVelocityToGrid() override;

    void applyBodyForces() override;

    void countParticles() override;

    void reseedParticles() override;

    void seedInitialFluid() override;

    void updateMaterials() override;

    void applyPressuresToVelocityField(std::vector<double> &pressures) override;

    void particleUpdate() override;

    void bumpParticles();

protected:

    double divergenceAt(int i, int j) override;

    float getFaceFractionUSample(int i, int j);
    float getFaceFractionVSample(int i, int j);

    double getWeightedDensityForUSample(int i, int j);
    double getWeightedDensityForVSample(int i, int j);

    double getWeightedVelocityUSample(int i, int j);
    double getWeightedVelocityVSample(int i, int j);

    StaggeredVelocityGrid m_airVelocityGrid;
    StaggeredVelocityGrid m_savedAirVelocityGrid;

    SdfGrid m_airSdf;
    Grid2d<int> m_airParticleCounts;
};

#endif // MULTIFLIPSOLVER_H
