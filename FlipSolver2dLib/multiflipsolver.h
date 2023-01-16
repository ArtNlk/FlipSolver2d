#ifndef MULTIFLIPSOLVER_H
#define MULTIFLIPSOLVER_H

#include "flipsolver2d.h"

class MultiflipSolver : public FlipSolver
{
public:
    MultiflipSolver(int extrapRadius = 1, bool vonNeumannNeighbors = false);

    void calcPressureRhs(std::vector<double> &rhs) override;

    DynamicUpperTriangularSparseMatrix getPressureProjectionMatrix() override;

    void step() override;

    void advect() override;

    void updateSdf() override;

    void combineSdf();

    void extrapolateSdf(Grid2d<float>& sdfGrid);

    void particleToGrid() override;

    void applyBodyForces() override;

    void countParticles() override;

    void reseedParticles() override;

    void seedInitialFluid() override;

    void updateMaterials() override;

    void applyPressuresToVelocityField(std::vector<double> &pressures) override;

    void particleUpdate() override;

    void bumpParticles();

protected:

    double getWeightedDensityForUSample(int i, int j);
    double getWeightedDensityForVSample(int i, int j);

    double getWeightedVelocityUSample(int i, int j);
    double getWeightedVelocityVSample(int i, int j);
};

#endif // MULTIFLIPSOLVER_H
