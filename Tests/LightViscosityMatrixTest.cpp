#include <cstdlib>

#include "catch2/catch_test_macros.hpp"
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "Eigen/Sparse"

#include "lightviscosityweights.h"
#include "materialgrid.h"
#include "staggeredvelocitygrid.h"
#include "viscositymodel.h"

using MatrixType = Eigen::SparseMatrix<double, Eigen::RowMajor>;

void populateMaterialGrid(MaterialGrid& grid)
{
    for(int i = 0; i < grid.sizeI(); i++)
    {
        for(int j = 0; j < grid.sizeJ(); j++)
        {
            int mat = std::rand() % 3;
            switch(mat)
            {
            case 0:
                grid.setAt(i,j,FluidMaterial::EMPTY);
                break;

            case 1:
                grid.setAt(i,j,FluidMaterial::FLUID);
                break;

            case 2:
                grid.setAt(i,j,FluidMaterial::SOLID);
                break;
            }
        }
    }
}

void populateFloatGrid(Grid2d<float>& grid)
{
    for(int i = 0; i < grid.sizeI(); i++)
    {
        for(int j = 0; j < grid.sizeJ(); j++)
        {
            grid.setAt(i,j, static_cast <float> (std::rand()) / (static_cast <float> (RAND_MAX/10)));
        }
    }
}

MatrixType getEigenMatrix(StaggeredVelocityGrid& velocityGrid,
          const Grid2d<float>& viscosityGrid,
          const MaterialGrid& materialGrid,
          const float dt,
          const float dx,
          const float density)
{
    const LinearIndexable2d& indexer = viscosityGrid;

    const size_t size = indexer.linearSize();

    Eigen::SparseMatrix<double,Eigen::RowMajor> output = Eigen::SparseMatrix<double>();
    output.resize(size,size);
    output.reserve(Eigen::VectorXi::Constant(size,10));

    const double scale = dt;
    //const double scale = 1.0;

    for(ssize_t i = 0; i < indexer.sizeI(); i++)
    {
        for(ssize_t j = 0; j < indexer.sizeJ(); j++)
        {
            ssize_t idx = indexer.linearIndex(i,j);

            double diag = 4.0;
            double ip1Neighbor = 1.0;
            double jp1Neighbor = 1.0;
            double im1Neighbor = 1.0;
            double jm1Neighbor = 1.0;

            if(materialGrid.isSolid(i,j))
            {
                output.coeffRef(idx,idx) = 1.0;
                continue;
            }

            ip1Neighbor *= (materialGrid.isSolid(i+1,j) ? viscosityGrid.at(i,j) : viscosityGrid.at(i+1,j)) * scale;
            jp1Neighbor *= (materialGrid.isSolid(i, j+1) ? viscosityGrid.at(i,j) : viscosityGrid.at(i,j+1)) * scale;
            im1Neighbor *= (materialGrid.isSolid(i-1, j) ? viscosityGrid.at(i,j) : viscosityGrid.at(i-1,j)) * scale;
            jm1Neighbor *= (materialGrid.isSolid(i, j-1) ? viscosityGrid.at(i,j) : viscosityGrid.at(i,j-1)) * scale;
            diag *= viscosityGrid.at(i,j) * scale;
            diag += 1.;

            output.coeffRef(idx,idx) = diag;
            if(indexer.inBounds(i+1,j))
            {
                output.coeffRef(idx,indexer.linearIndex(i+1,j)) = ip1Neighbor;
                output.coeffRef(indexer.linearIndex(i+1,j),idx) = ip1Neighbor;
            }

            if(indexer.inBounds(i,j+1))
            {
                output.coeffRef(idx,indexer.linearIndex(i,j+1)) = jp1Neighbor;
                output.coeffRef(indexer.linearIndex(i,j+1),idx) = jp1Neighbor;
            }

            if(indexer.inBounds(i-1,j))
            {
                output.coeffRef(idx,indexer.linearIndex(i-1,j)) = im1Neighbor;
                output.coeffRef(indexer.linearIndex(i-1,j),idx) = im1Neighbor;
            }

            if(indexer.inBounds(i,j-1))
            {
                output.coeffRef(idx,indexer.linearIndex(i,j-1)) = jm1Neighbor;
                output.coeffRef(indexer.linearIndex(i,j-1),idx) = jm1Neighbor;
            }
        }
    }

    output.makeCompressed();

    return output;
}

void populateInputVectors(std::vector<double>& myVec, Eigen::VectorXd& eigenVec)
{
    for(int i = 0; i < myVec.size(); i++)
    {
        double val = std::rand();
        myVec[i] = val;
        eigenVec.coeffRef(i) = val;
    }
};

TEST_CASE("Viscosity matrix matches Eigen")
{
    const int sizeI = 256;
    const int sizeJ = 256;

    StaggeredVelocityGrid dummy(sizeI, sizeJ);
    Grid2d<float> viscosityGrid(sizeI, sizeJ);
    MaterialGrid materials(sizeI, sizeJ, FluidMaterial::EMPTY);

    float dt = 0.1f;
    float dx = 1.f/256.f;
    float density = 1.f;

    LightViscosityWeights weights = LightViscosityModel::getMatrix(dummy, viscosityGrid, materials, dt, dx, density);
    MatrixType eigenMat = getEigenMatrix(dummy, viscosityGrid, materials, dt, dx, density);

    std::vector<double> vIn(materials.linearSize());
    Eigen::VectorXd eigenVin;

    std::vector<double> vOut(materials.linearSize());
    Eigen::VectorXd eigenVout;

    eigenVin.resize(vIn.size());
    eigenVout.resize(vOut.size());

    populateInputVectors(vIn, eigenVin);

    weights.multiply(vIn, vOut);

    eigenVout = eigenMat * eigenVin;

    for(int i = 0; i < vOut.size(); i++)
    {
        CHECK_THAT(vOut[i], Catch::Matchers::WithinAbs(eigenVout.coeff(i), 0.0000001));
    }
}
