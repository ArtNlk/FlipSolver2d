#include "catch2/catch_test_macros.hpp"
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <random>
#include <vector>

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

    for(size_t i = 0; i < materialGrid.sizeI(); i++)
    {
        for(size_t j = 0; j < materialGrid.sizeJ(); j++)
        {
            const size_t linIdx = indexer.linearIndex(i,j);

            if(!materialGrid.isFluid(i,j))
            {
                if(linIdx >= threadRanges.at(currRangeIdx).end)
                {
                    output.endThreadDataRange();
                    currRangeIdx++;
                }
                continue;
            }

            IndexedIPPCoefficientUnit unit;
            unit.unitIndex = linIdx;

            const ssize_t iNegLinIdx = indexer.linearIdxOfOffset(linIdx,-1,0);
            const ssize_t iPosLinIdx = indexer.linearIdxOfOffset(linIdx,1,0);
            const ssize_t jNegLinIdx = indexer.linearIdxOfOffset(linIdx,0,-1);
            const ssize_t jPosLinIdx = indexer.linearIdxOfOffset(linIdx,0,1);

            unit.iNeg = 1.0/(materialGrid.nonsolidNeighborCount(iNegLinIdx)*scale);
            unit.iPos = 1.0/(materialGrid.nonsolidNeighborCount(iPosLinIdx)*scale);
            unit.jNeg = 1.0/(materialGrid.nonsolidNeighborCount(jNegLinIdx)*scale);
            unit.jPos = 1.0/(materialGrid.nonsolidNeighborCount(jPosLinIdx)*scale);

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
    const int sizeI = 256;
    const int sizeJ = 256;
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

    for(int i = 0; i < myOutput.size(); i++)
    {
        CHECK_THAT(myOutput[i], Catch::Matchers::WithinAbs(eigenOutput.coeff(i), 0.0000001));
    }
}
