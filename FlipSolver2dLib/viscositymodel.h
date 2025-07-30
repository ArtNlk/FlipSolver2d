#ifndef VISCOSITYMODEL_H
#define VISCOSITYMODEL_H

#include "grid2d.h"
#include "linearsolver.h"
#include "materialgrid.h"
#include "staggeredvelocitygrid.h"
#include "lightviscosityweights.h"
#include "heavyviscosityweights.h"

#include <Eigen/Sparse>

class ViscosityModel
{
public:
    using MatrixType = Eigen::SparseMatrix<double, Eigen::RowMajor>;

    ViscosityModel() = default;

    virtual ~ViscosityModel() = default;

    virtual int apply(StaggeredVelocityGrid& velocityGrid,
                     const Grid2d<float>& viscosityGrid,
                     const MaterialGrid& materialGrid,
                     float dt,
                     float dx,
                     float density) = 0;
protected:
    LinearSolver m_solver;
};

class LightViscosityModel : public ViscosityModel
{
public:
    int apply(StaggeredVelocityGrid &velocityGrid,
              const Grid2d<float> &viscosityGrid,
              const MaterialGrid &materialGrid,
              float dt,
              float dx,
              float density) override;

    static LightViscosityWeights getMatrix(StaggeredVelocityGrid& velocityGrid,
                         const Grid2d<float>& viscosityGrid,
                         const MaterialGrid& materialGrid,
                         const float dt,
                         const float dx,
                         const float density);

protected:
    void fillRhs(std::vector<double>& rhs,
                 const Grid2d<float>& velocityGrid,
                 const LinearIndexable2d& indexer,
                 float density);

    void applyResult(Grid2d<float>& velocityGrid, const LinearIndexable2d& indexer, const std::vector<double>& result, float density);
};

class HeavyViscosityModel : public ViscosityModel
{
    int apply(StaggeredVelocityGrid& velocityGrid,
               const Grid2d<float>& viscosityGrid,
               const MaterialGrid& materialGrid,
               float dt,
               float dx,
               float density) override;

    HeavyViscosityWeights getMatrix(StaggeredVelocityGrid& velocityGrid,
                   const Grid2d<float>& viscosityGrid,
                   const MaterialGrid& materialGrid,
                   float dt,
                   float dx,
                   float density);

    void fillRhs(std::vector<double>& rhs, const StaggeredVelocityGrid& velocityGrid, float density);

    void applyResult(StaggeredVelocityGrid& velocityGrid, const std::vector<double>& result);
};
#endif // VISCOSITYMODEL_H
