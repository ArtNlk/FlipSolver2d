#include "catch2/catch_test_macros.hpp"
#include "catch2/matchers/catch_matchers_floating_point.hpp"

#include <cstdlib>
#include <vector>
#include "Eigen/Sparse"

#include "staticmatrix.h"

const size_t matRowSize = 9;
const size_t matrixSize = 1024;

void populateMatricies(DynamicMatrix<matRowSize>& myMatrix, Eigen::SparseMatrix<double,Eigen::RowMajor>& eigenMatrix)
{
    eigenMatrix.reserve(matrixSize * matRowSize);

    for(size_t rowIndex = 0; rowIndex < matrixSize; rowIndex++)
    {
        const float rowProb = static_cast <float> (std::rand()) / static_cast <float> (RAND_MAX);
        if(rowProb > 0.7)
        {
            continue;
        }

        for(size_t colIndex = 0; colIndex < matRowSize; colIndex++)
        {
            const float colProb = static_cast <float> (std::rand()) / static_cast <float> (RAND_MAX);
            if(colProb > 0.7)
            {
                continue;
            }

            const double value = static_cast <double> (std::rand()) / (static_cast <double> (RAND_MAX/10));

            myMatrix.addValue(rowIndex, colIndex, value);
            eigenMatrix.coeffRef(rowIndex, colIndex) = value;
        }
    }
}

static void populateVectors(std::vector<double>& myVector, Eigen::VectorXd& eigenVector)
{
    myVector.resize(matrixSize, 0.0);
    eigenVector.resize(matrixSize);

    for(size_t idx = 0; idx < matrixSize; idx++)
    {
        const double value = static_cast <double> (std::rand()) / (static_cast <double> (RAND_MAX/10));

        myVector[idx] = value;
        eigenVector[idx] = value;
    }
}

TEST_CASE("Matrix impl matches Eigen") {
    DynamicMatrix<matRowSize> myMatrixSource(matrixSize);
    Eigen::SparseMatrix<double,Eigen::RowMajor> eigenMatrix(matrixSize, matrixSize);
    populateMatricies(myMatrixSource, eigenMatrix);

    StaticMatrix myMatrix = StaticMatrix::fromDynamic(myMatrixSource);

    std::vector<double> myVector;
    Eigen::VectorXd eigenVector;
    populateVectors(myVector, eigenVector);

    std::vector<double> myResult(matrixSize, 0.0);
    Eigen::VectorXd eigenResult(matrixSize);
    eigenResult.fill(0.0);

    myMatrix.multiply(myVector, myResult);
    eigenResult = eigenMatrix*eigenVector;

    for(size_t idx = 0; idx < matrixSize; idx++)
    {
        CHECK_THAT(myResult[idx], Catch::Matchers::WithinAbs(eigenResult.coeff(idx), 0.0000001));
    }
}
