#include "catch2/catch_test_macros.hpp"
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <random>
#include <vector>
#include <array>

#include "InversePoissonPreconditioner.h"
#include "materialgrid.h"

static std::random_device rd;
static std::mt19937 gen(rd());

InversePoissonPreconditioner getCustomPrecond(double stepDt, double density, double dx, const MaterialGrid& materialGrid)
{
    const double scale = stepDt / (density * dx * dx);

    InversePoissonPreconditioner output(materialGrid.linearSize()*0.33, materialGrid);

    const LinearIndexable2d& indexer = static_cast<LinearIndexable2d>(materialGrid);

    std::vector<Range> threadRanges = ThreadPool::i()->splitRange(materialGrid.linearSize());
    size_t currRangeIdx = 0;

    std::vector<std::array<double,2>> tempData(materialGrid.linearSize(),{0.0,0.0});

    for(size_t i = 0; i < materialGrid.sizeI(); i++)
    {
        for(size_t j = 0; j < materialGrid.sizeJ(); j++)
        {
            const size_t linIdx = indexer.linearIndex(i,j);

            if(!materialGrid.isFluid(i,j))
            {
                tempData[linIdx] = {0.0,0.0};
            }

            double diag = materialGrid.nonsolidNeighborCount(i,j);
            double iNeg = materialGrid.isFluid(i-1, j) ? scale : 0.0;
            double jNeg = materialGrid.isFluid(i, j-1) ? scale : 0.0;

            tempData[linIdx] = {1.0-(jNeg/diag), 1.0-(iNeg/diag)};
        }
    }

    for(size_t i = 0; i < materialGrid.sizeI(); i++)
    {
        for(size_t j = 0; j < materialGrid.sizeJ(); j++)
        {
            const ssize_t linIdx = indexer.linearIndex(i,j);

            if(!materialGrid.isFluid(i,j))
            {
                if(linIdx >= threadRanges.at(currRangeIdx).end)
                {
                    output.endThreadDataRange();
                    currRangeIdx++;
                }
                continue;
            }

            std::array<double,2> currRowData= tempData[linIdx];

            IndexedIPPCoefficientUnit unit;
            unit.unitIndex = linIdx;

            ssize_t b0Idx = linIdx - indexer.iLinearOffset();
            ssize_t b1Idx = linIdx - indexer.iLinearOffset() + indexer.jLinearOffset();
            ssize_t b2Idx = linIdx - indexer.jLinearOffset();
            ssize_t b3Idx = linIdx;
            ssize_t b4Idx = linIdx + indexer.jLinearOffset();
            ssize_t b5Idx = linIdx + indexer.iLinearOffset() - indexer.jLinearOffset();
            ssize_t b6Idx = linIdx + indexer.iLinearOffset();

            double b1data = indexer.inBounds(b1Idx) ? tempData[b1Idx][1] : 1.0;
            double b4data = indexer.inBounds(b4Idx) ? tempData[b4Idx][1] : 1.0;
            double b5data = indexer.inBounds(b5Idx) ? tempData[b5Idx][0] : 1.0;
            double b6data = indexer.inBounds(b6Idx) ? tempData[b6Idx][0] : 1.0;

            unit.data[0] = currRowData[0];
            unit.data[1] = currRowData[0] * b1data;
            unit.data[2] = currRowData[1];
            unit.data[3] = currRowData[0] * currRowData[0] + currRowData[1]*currRowData[1];
            unit.data[4] = b4data;
            unit.data[5] = currRowData[1] * b5data;
            unit.data[6] = b6data;

            if(linIdx >= threadRanges.at(currRangeIdx).end)
            {
                output.endThreadDataRange();
                currRangeIdx++;
            }

            output.add(unit);
        }
    }

    if(currRangeIdx != threadRanges.size())
    {
        output.endThreadDataRange();
    }

    return output;
}

Eigen::SparseMatrix<double,Eigen::RowMajor> getEigenPressureProjectionMatrix(double stepDt,
                                                                              double fluidDensity,
                                                                              double dx,
                                                                              const MaterialGrid& materialGrid)
{
    Eigen::SparseMatrix<double,Eigen::RowMajor> output = Eigen::SparseMatrix<double>();
    output.resize(materialGrid.linearSize(),materialGrid.linearSize());
    output.reserve(Eigen::VectorXi::Constant(materialGrid.linearSize(),6));

    const double scale = stepDt / (fluidDensity * dx * dx);

    const LinearIndexable2d& indexer = materialGrid;

    for(int i = 0; i < materialGrid.sizeI(); i++)
    {
        for(int j = 0; j < materialGrid.sizeJ(); j++)
        {
            const int linIdx = indexer.linearIndex(i,j);
            if(!materialGrid.isFluid(i,j))
            {
                output.insert(linIdx,linIdx) = 1.f;
                continue;
            }

            const int linIdxAx = indexer.linearIdxOfOffset(linIdx,1,0);
            const int linIdxAy = indexer.linearIdxOfOffset(linIdx,0,1);

            double diag = 0.0;
            //X Neighbors
            if(materialGrid.isFluid(i-1,j))
            {
                diag += scale;
            }else if(materialGrid.isEmpty(i-1,j))
            {
                diag += scale;
            }

            if(materialGrid.isFluid(i+1,j))
            {
                diag += scale;
                if(materialGrid.inBounds(linIdxAx))
                {
                    output.insert(linIdxAx,linIdx) = -scale;
                    output.insert(linIdx,linIdxAx) = -scale;
                }
            } else if(materialGrid.isEmpty(i+1,j))
            {
                diag += scale;
            }

            //Y Neighbors
            if(materialGrid.isFluid(i,j-1))
            {
                diag += scale;
            }else if(materialGrid.isEmpty(i,j-1))
            {
                diag += scale;
            }

            if(materialGrid.isFluid(i,j+1))
            {
                diag += scale;
                if(materialGrid.inBounds(linIdxAy))
                {
                    output.insert(linIdx,linIdxAy) = -scale;
                    output.insert(linIdxAy,linIdx) = -scale;
                }
            } else if(materialGrid.isEmpty(i,j+1))
            {
                diag += scale;
            }

            output.insert(linIdx,linIdx) = diag;
        }
    }

    //    for(int i = 0; i <  m_sizeI; i++)
    //    {
    //        for(int j = 0; j <  m_sizeJ; j++)
    //        {
    //            if(m_materialGrid.isFluid(i,j))
    //            {
    //                //X Neighbors
    //                if(m_materialGrid.isFluid(i-1,j))
    //                {
    //                    output.addToAdiag(i,j,scale,indexer);
    //                }else if(m_materialGrid.isEmpty(i-1,j))
    //                {
    //                    output.addToAdiag(i,j,scale,  indexer);
    //                }

    //                if(m_materialGrid.isFluid(i+1,j))
    //                {
    //                    output.addToAdiag(i,j,scale,  indexer);
    //                    output.setAx(i,j,-scale,  indexer);
    //                } else if(m_materialGrid.isEmpty(i+1,j))
    //                {
    //                    output.addToAdiag(i,j,scale,  indexer);
    //                }

    //                //Y Neighbors
    //                if(m_materialGrid.isFluid(i,j-1))
    //                {
    //                    output.addToAdiag(i,j,scale,  indexer);
    //                }else if(m_materialGrid.isEmpty(i,j-1))
    //                {
    //                    output.addToAdiag(i,j,scale,  indexer);
    //                }

    //                if(m_materialGrid.isFluid(i,j+1))
    //                {
    //                    output.addToAdiag(i,j,scale,  indexer);
    //                    output.setAy(i,j,-scale,  indexer);
    //                } else if(m_materialGrid.isEmpty(i,j+1))
    //                {
    //                    output.addToAdiag(i,j,scale,  indexer);
    //                }
    //            }
    //        }
    //    }

    output.makeCompressed();

    return output;
}

Eigen::SparseMatrix<double,Eigen::RowMajor> factorize(const Eigen::SparseMatrix<double,Eigen::RowMajor>& mat)
{
    Eigen::SparseMatrix<double,Eigen::RowMajor> output;
    auto identity = mat;
    identity.setIdentity();
    Eigen::VectorXd invdiag;
    invdiag.resize(mat.cols());

    for(ssize_t j=0; j<mat.outerSize(); ++j)
    {
        Eigen::SparseMatrix<double,Eigen::RowMajor>::InnerIterator it(mat,j);
        while(it && it.index()!=j)
        {
            ++it;
        }

        if(it && it.index()==j && it.value()!=double(0))
        {
            invdiag(j) = double(1)/it.value();
        }
        else
        {
            invdiag(j) = double(1);
        }
    }
    output = mat.template triangularView<Eigen::StrictlyLower>();

    for(ssize_t j=0; j<output.outerSize(); ++j)
    {
        Eigen::SparseMatrix<double,Eigen::RowMajor>::InnerIterator it(output,j);
        for(;it && it.index() < j; ++it)
        {
            it.valueRef() *= invdiag(it.col());
        }
    }
    output = identity - output;
    output = output * output.transpose();

    return output;
}

MaterialGrid getMaterialGrid(int sizeI, int sizeJ)
{
    std::uniform_int_distribution<> distr(0, 2);

    MaterialGrid output(sizeI, sizeJ, FluidMaterial::EMPTY);

    for(int i = 0; i < output.sizeI(); i++)
    {
        for(int j = 0; j < output.sizeJ(); j++)
        {
            FluidMaterial m;
            switch(distr(gen))
            {
            case 0:
                m = FluidMaterial::EMPTY;
                break;
            case 1:
                m = FluidMaterial::SOLID;
                break;
            case 2:
            default:
                m = FluidMaterial::FLUID;
                break;
            }

            output.setAt(i,j, m);
        }
    }

    return output;
}

void populateVectors(std::vector<double>& myVec, Eigen::VectorXd& eigenVec)
{
    std::uniform_real_distribution dist(0.0,100000.0);

    for(int i = 0; i < myVec.size(); i++)
    {
        double val = dist(gen);
        myVec.at(i) = val;
        eigenVec.coeffRef(i) = val;
    }
}

TEST_CASE("Incomplete poisson preconditioner matches Eigen")
{
    const int sizeI = 64;
    const int sizeJ = 64;
    const double stepDt = 0.03;
    const double density = 0.1;
    const double dx = 0.1;

    const int linearSize = sizeI*sizeJ;

    MaterialGrid sourceGrid = getMaterialGrid(sizeI, sizeJ);

    InversePoissonPreconditioner precond = getCustomPrecond(stepDt, density, dx, sourceGrid);
    Eigen::SparseMatrix<double,Eigen::RowMajor> sourceMat = getEigenPressureProjectionMatrix(stepDt, density, dx, sourceGrid);
    Eigen::SparseMatrix<double,Eigen::RowMajor> eigenPrecond = factorize(sourceMat);

    std::vector<double> myVec(linearSize);

    Eigen::VectorXd eigenVec;
    eigenVec.resize(linearSize);

    populateVectors(myVec, eigenVec);

    std::vector<double> myOutput(linearSize, 0.0);

    Eigen::VectorXd eigenOutput;
    eigenOutput.resize(linearSize);
    eigenOutput.fill(0.0);

    precond.multiply(myVec, myOutput);

    eigenOutput = eigenPrecond * eigenVec;

    // std::cout << "============MAT_START=============" << std::endl;
    // for(int k = 0; k < eigenPrecond.outerSize(); ++k) {
    //     for(Eigen::SparseMatrix<double,Eigen::RowMajor>::InnerIterator it(eigenPrecond,k);it;++it)
    //     {
    //         std::cout << it.row() << ' ' << it.col() << ' ' << it.value() << '\n';
    //     }
    // }

    // std::cout << "=============MAT_END==============" << std::endl;

    for(int i = 0; i < myOutput.size(); i++)
    {
        CHECK_THAT(myOutput[i], Catch::Matchers::WithinAbs(eigenOutput.coeff(i), 0.0000001));
    }
}
