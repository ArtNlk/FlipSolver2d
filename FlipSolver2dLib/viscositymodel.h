#ifndef VISCOSITYMODEL_H
#define VISCOSITYMODEL_H

#include "grid2d.h"
#include "linearindexable2d.h"
#include "materialgrid.h"
#include <Eigen/Sparse>

class ViscosityModel
{
public:
    using MatrixType = Eigen::SparseMatrix<double, Eigen::RowMajor>;

    ViscosityModel() = default;

    virtual ~ViscosityModel() = default;

    virtual MatrixType getMatrix(const LinearIndexable2d& indexerU,
                                 const LinearIndexable2d& indexerV,
                                 const Grid2d<float>& viscosityGrid,
                                 const MaterialGrid& materialGrid,
                                 const float dt,
                                 const float dx,
                                 const float density) = 0;
};

class LightViscosityModel : public ViscosityModel
{
    MatrixType getMatrix(const LinearIndexable2d& indexerU,
                        const LinearIndexable2d& indexerV,
                        const Grid2d<float>& viscosityGrid,
                        const MaterialGrid& materialGrid,
                        const float dt,
                        const float dx,
                        const float density) override;
};

class HeavyViscosityModel : public ViscosityModel
{
    MatrixType getMatrix(const LinearIndexable2d& indexerU,
                         const LinearIndexable2d& indexerV,
                         const Grid2d<float>& viscosityGrid,
                         const MaterialGrid& materialGrid,
                         const float dt,
                         const float dx,
                         const float density) override;
};
#endif // VISCOSITYMODEL_H
