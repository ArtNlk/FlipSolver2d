#ifndef VISCOSITYMODEL_H
#define VISCOSITYMODEL_H

#include "grid2d.h"
#include "materialgrid.h"
#include "staggeredvelocitygrid.h"
#include <Eigen/Sparse>

class ViscosityModel
{
public:
    using MatrixType = Eigen::SparseMatrix<double, Eigen::RowMajor>;

    ViscosityModel() = default;

    virtual ~ViscosityModel() = default;

    virtual void apply(StaggeredVelocityGrid& velocityGrid,
                     const Grid2d<float>& viscosityGrid,
                     const MaterialGrid& materialGrid,
                     float dt,
                     float dx,
                     float density) = 0;
protected:
    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Upper> m_viscositySolver;
};

class LightViscosityModel : public ViscosityModel
{
    void apply(StaggeredVelocityGrid& velocityGrid,
               const Grid2d<float>& viscosityGrid,
               const MaterialGrid& materialGrid,
               float dt,
               float dx,
               float density) override;

    MatrixType getMatrix(StaggeredVelocityGrid& velocityGrid,
                         const Grid2d<float>& viscosityGrid,
                         const MaterialGrid& materialGrid,
                         const float dt,
                         const float dx,
                         const float density);

    void fillRhs(Eigen::VectorXd& rhs,
                 const Grid2d<float>& velocityGrid,
                 const LinearIndexable2d& indexer,
                 float density);

    void applyResult(Grid2d<float>& velocityGrid, const LinearIndexable2d& indexer, const Eigen::VectorXd& result, float density);
};

class HeavyViscosityModel : public ViscosityModel
{
    void apply(StaggeredVelocityGrid& velocityGrid,
               const Grid2d<float>& viscosityGrid,
               const MaterialGrid& materialGrid,
               float dt,
               float dx,
               float density) override;

    MatrixType getMatrix(StaggeredVelocityGrid& velocityGrid,
                   const Grid2d<float>& viscosityGrid,
                   const MaterialGrid& materialGrid,
                   float dt,
                   float dx,
                   float density);

    Eigen::VectorXd getRhs(const StaggeredVelocityGrid& velocityGrid, float density);

    void applyResult(StaggeredVelocityGrid& velocityGrid, const Eigen::VectorXd& result);
};
#endif // VISCOSITYMODEL_H
