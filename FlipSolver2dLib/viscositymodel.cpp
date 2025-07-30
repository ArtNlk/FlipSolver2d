#include "viscositymodel.h"
#include "Eigen/src/Core/Matrix.h"
#include "IdentityPreconditioner.h"
#include "lightviscosityweights.h"
#include "linearindexable2d.h"

int LightViscosityModel::apply(StaggeredVelocityGrid& velocityGrid,
                              const Grid2d<float>& viscosityGrid,
                              const MaterialGrid& materialGrid,
                              float dt,
                              float dx,
                              float density)
{
    std::vector<double> rhs;
    std::vector<double> result;
    rhs.resize(viscosityGrid.linearSize());

    result.resize(viscosityGrid.linearSize());

    auto viscosityMatrix = getMatrix(velocityGrid,
                                     viscosityGrid,
                                     materialGrid,
                                     dt,
                                     dx,
                                     density);

    fillRhs(rhs,velocityGrid.velocityGridU(),viscosityGrid,density);

    IdentityPreconditioner precond;

    int iters = 0;
    iters = m_solver.solve(viscosityMatrix,precond,result,rhs,600,1e-6);
    if(iters == 600) {
        std::cout << "Viscosity solver U solving failed!\n";
        return -1;
    }
    std::cout << "Viscosity U done with " << iters << " iterations\n" << std::endl;

    applyResult(velocityGrid.velocityGridU(), viscosityGrid, result, density);

    fillRhs(rhs,velocityGrid.velocityGridV(),viscosityGrid,density);
    iters = m_solver.solve(viscosityMatrix,precond,result,rhs,600,1e-6);
    if(iters == 600) {
        std::cout << "Viscosity solver V solving failed!\n";
        return -1;
    }
    std::cout << "Viscosity V done with " << iters << " iterations\n";

    applyResult(velocityGrid.velocityGridV(), viscosityGrid, result, density);

    return iters;

    // if(anyNanInf(velocityGrid.velocityGridU().data()))
    // {
    //     std::cout << "NaN or inf in U velocity after viscosity!\n" << std::flush;
    // }

    // if(anyNanInf(velocityGrid.velocityGridV().data()))
    // {
    //     std::cout << "NaN or inf in V velocity after viscosity!\n" << std::flush;
    // }
}

LightViscosityWeights LightViscosityModel::getMatrix(StaggeredVelocityGrid &velocityGrid,
                                                          const Grid2d<float> &viscosityGrid,
                                                          const MaterialGrid &materialGrid,
                                                          float dt,
                                                          float dx,
                                                          float density)
{
    const LinearIndexable2d& indexer = viscosityGrid;

    const size_t size = indexer.linearSize();

    LightViscosityWeights output(materialGrid);

    const double scale = dt;
    //const double scale = 1.0;

    std::vector<Range> threadRanges = ThreadPool::i()->splitRange(size);
    size_t currRangeIdx = 0;

    for(ssize_t i = 0; i < indexer.sizeI(); i++)
    {
        for(ssize_t j = 0; j < indexer.sizeJ(); j++)
        {
            ssize_t idx = indexer.linearIndex(i,j);

            LightViscosityWeightsUnit unit;

            double diag = 4.0;
            double ip1Neighbor = 1.0;
            double jp1Neighbor = 1.0;
            double im1Neighbor = 1.0;
            double jm1Neighbor = 1.0;

            if(materialGrid.isSolid(i,j))
            {
                unit.diag = 1.0;
                output.add(unit);
                continue;
            }

            im1Neighbor *= (materialGrid.isSolid(i-1, j) ? viscosityGrid.at(i,j) : viscosityGrid.at(i-1,j)) * scale;
            jm1Neighbor *= (materialGrid.isSolid(i, j-1) ? viscosityGrid.at(i,j) : viscosityGrid.at(i,j-1)) * scale;
            ip1Neighbor *= (materialGrid.isSolid(i+1,j) ? viscosityGrid.at(i,j) : viscosityGrid.at(i+1,j)) * scale;
            jp1Neighbor *= (materialGrid.isSolid(i, j+1) ? viscosityGrid.at(i,j) : viscosityGrid.at(i,j+1)) * scale;

            diag *= viscosityGrid.at(i,j) * scale;
            diag += 1.;

            unit.diag = diag;

            if(materialGrid.inBounds(i-1,j))
            {
                unit.im1 = im1Neighbor;
            }

            if(materialGrid.inBounds(i,j-1))
            {
                unit.jm1 = jm1Neighbor;
            }

            if(materialGrid.inBounds(i+1,j))
            {
                unit.im1 = ip1Neighbor;
            }

            if(materialGrid.inBounds(i,j+1))
            {
                unit.im1 = jp1Neighbor;
            }

            if(idx >= threadRanges.at(currRangeIdx).end)
            {
                output.endThreadDataRange();
                currRangeIdx++;
            }

            output.add(unit);
        }
    }

    if(currRangeIdx < threadRanges.size())
    {
        output.endThreadDataRange();
    }

    return output;
}

void LightViscosityModel::fillRhs(std::vector<double>& rhs, const Grid2d<float> &velocityGrid, const LinearIndexable2d &indexer, float density)
{
    for (ssize_t i = 0; i < indexer.sizeI(); i++)
    {
        for (ssize_t j = 0; j < indexer.sizeJ(); j++)
        {
            ssize_t idx = indexer.linearIndex(i,j);
            rhs[idx] = density * velocityGrid.at(i,j);
        }
    }
}

void LightViscosityModel::applyResult(Grid2d<float> &velocityGrid, const LinearIndexable2d &indexer, const std::vector<double> &result, float density)
{
    for (ssize_t i = 0; i < indexer.sizeI(); i++)
    {
        for (ssize_t j = 0; j < indexer.sizeJ(); j++)
        {
            ssize_t idx = indexer.linearIndex(i,j);
            velocityGrid.setAt(i,j,result[idx]/density);
        }
    }
}

int HeavyViscosityModel::apply(StaggeredVelocityGrid& velocityGrid,
                              const Grid2d<float>& viscosityGrid,
                              const MaterialGrid& materialGrid,
                              float dt,
                              float dx,
                              float density)
{
    std::vector<double> rhs;
    fillRhs(rhs,velocityGrid,density);

    std::vector<double> result;
    result.resize(rhs.size());

    auto viscosityMatrix = getMatrix(velocityGrid,
                                     viscosityGrid,
                                     materialGrid,
                                     dt,
                                     dx,
                                     density);

    IdentityPreconditioner precond;

    int iters = m_solver.solve(viscosityMatrix,precond,result,rhs,600,1e-6);

    if(iters == 600) {
        std::cout << "Viscosity solver combined solving failed!\n";
        return -1;
    }
    std::cout << "Viscosity combined done with " << iters << " iterations\n";

    applyResult(velocityGrid, result);

    return iters;
    // if(anyNanInf(velocityGrid.velocityGridU().data()))
    // {
    //     std::cout << "NaN or inf in U velocity after viscosity!\n" << std::flush;
    // }

    // if(anyNanInf(velocityGrid.velocityGridV().data()))
    // {
    //     std::cout << "NaN or inf in V velocity after viscosity!\n" << std::flush;
    // }
}

HeavyViscosityWeights HeavyViscosityModel::getMatrix(StaggeredVelocityGrid &velocityGrid,
                                    const Grid2d<float> &viscosityGrid,
                                    const MaterialGrid &materialGrid,
                                    float dt,
                                    float dx,
                                    float density)
{
    const LinearIndexable2d& indexer = viscosityGrid;
    const LinearIndexable2d& indexerU = velocityGrid.velocityGridU();
    const LinearIndexable2d& indexerV = velocityGrid.velocityGridV();

    const size_t size = indexerU.linearSize() +
                        indexerV.linearSize();

    HeavyViscosityWeights output(indexer.linearSize() * 2 * 0.5);

    const float scaleTwoDt = 2*dt / (dx * dx);
    const float scaleTwoDx = dt / (2 * dx * dx);
    const ssize_t vBaseIndex = indexerU.linearSize();

    for(ssize_t i = 0; i < indexer.sizeI()+1; i++)
    {
        for(ssize_t j = 0; j < indexer.sizeJ()+1; j++)
        {
            ssize_t idxU = indexerU.linearIndex(i,j);
            ssize_t idxV = indexerV.linearIndex(i,j);

            if(idxU != -1)
            {
                float fi = static_cast<float>(i);
                float fj = static_cast<float>(j);

                //U component
                output.data()[idxU]idxU) += density;

                ssize_t uImOneLinearIdx = indexerU.linearIndex(i-1,j);

                if(uImOneLinearIdx != -1)
                {
                    output.data()[idxU]
                                    uImOneLinearIdx) += -scaleTwoDt * viscosityGrid.getAt(i-1,j);

                    output.data()[idxU]
                                    idxU) += scaleTwoDt * viscosityGrid.getAt(i-1,j);
                }

                ssize_t uIpOneLinearIdx = indexerU.linearIndex(i+1,j);

                if(uIpOneLinearIdx != -1)
                {
                    output.data()[idxU]
                                    uIpOneLinearIdx) += -scaleTwoDt * viscosityGrid.getAt(i,j);

                    output.data()[idxU]
                                    idxU) += scaleTwoDt * viscosityGrid.getAt(i,j);
                }

                ssize_t uJmOneLinearIdx = indexerU.linearIndex(i,j-1);
                ssize_t vImOneLinearIdx = indexerV.linearIndex(i-1,j);

                if(uJmOneLinearIdx != -1
                    && idxV != -1
                    && vImOneLinearIdx != -1)
                {
                    float lerpedViscosity = viscosityGrid.interpolateAt(fi-0.5f,fj-0.5f);
                    //lerpedViscosity = tempVisc;
                    output.data()[idxU]
                                    uJmOneLinearIdx) += -scaleTwoDx * lerpedViscosity;

                    output.data()[idxU]
                                    vBaseIndex + idxV) += scaleTwoDx * lerpedViscosity;

                    output.data()[idxU]
                                    vBaseIndex + vImOneLinearIdx) += -scaleTwoDx * lerpedViscosity;

                    output.data()[idxU]idxU) += scaleTwoDx * lerpedViscosity;
                }

                ssize_t uJpOneLinearIdx = indexerU.linearIndex(i,j+1);
                ssize_t vJpOneLinearIdx = indexerV.linearIndex(i,j+1);
                ssize_t vImOneJpOneLinearIdx = indexerV.linearIndex(i-1,j+1);

                if(uJpOneLinearIdx != -1
                    && vJpOneLinearIdx != -1
                    && vImOneJpOneLinearIdx != -1)
                {
                    float lerpedViscosity = viscosityGrid.interpolateAt(fi-0.5f,fj+0.5f);
                    //lerpedViscosity = tempVisc;
                    output.data()[idxU]
                                    uJpOneLinearIdx) += -scaleTwoDx * lerpedViscosity;

                    output.data()[idxU]
                                    vBaseIndex + vJpOneLinearIdx) += -scaleTwoDx * lerpedViscosity;

                    output.data()[idxU]
                                    vBaseIndex + vImOneJpOneLinearIdx) += scaleTwoDx * lerpedViscosity;

                    output.data()[idxU]idxU) += scaleTwoDx * lerpedViscosity;
                }
            }

            //V component
            if(idxV != -1)
            {
                float fi = static_cast<float>(i);
                float fj = static_cast<float>(j);
                ssize_t vBaseIndex = indexerU.linearSize();
                output.data()[vBaseIndex + idxV]vBaseIndex + idxV) += density;

                ssize_t vJmOneLinearIdx = indexerV.linearIndex(i,j-1);

                if(vJmOneLinearIdx != -1)
                {
                    output.data()[vBaseIndex + idxV]
                                    vBaseIndex + vJmOneLinearIdx) += -scaleTwoDt * viscosityGrid.getAt(i,j-1);

                    output.data()[vBaseIndex + idxV]
                                    vBaseIndex + idxV) += scaleTwoDt * viscosityGrid.getAt(i,j-1);
                }

                ssize_t vJpOneLinearIdx = indexerV.linearIndex(i,j+1);

                if(vJpOneLinearIdx != -1)
                {
                    output.data()[vBaseIndex + idxV]
                                    vBaseIndex + vJpOneLinearIdx)
                        += -scaleTwoDt * viscosityGrid.getAt(i,j);

                    output.data()[vBaseIndex + idxV]
                                    vBaseIndex + idxV)
                        += scaleTwoDt * viscosityGrid.getAt(i,j);
                }

                ssize_t uJmOneLinearIdx = indexerU.linearIndex(i,j-1);
                ssize_t vImOneLinearIdx = indexerV.linearIndex(i-1,j);

                if(idxU != -1
                    && uJmOneLinearIdx != -1
                    && vImOneLinearIdx != -1)
                {
                    float lerpedViscosity = viscosityGrid.interpolateAt(fi-0.5f,fj-0.5f);
                    //lerpedViscosity = tempVisc;

                    output.data()[vBaseIndex + idxV]
                                    idxU) +=
                        scaleTwoDx * lerpedViscosity;

                    output.data()[vBaseIndex + idxV]
                                    uJmOneLinearIdx) +=
                        -scaleTwoDx * lerpedViscosity;

                    output.data()[vBaseIndex + idxV]
                                    vBaseIndex + vImOneLinearIdx) +=
                        -scaleTwoDx * lerpedViscosity;

                    output.data()[vBaseIndex + idxV]
                                    vBaseIndex + idxV) +=
                        scaleTwoDx * lerpedViscosity;
                }

                ssize_t uIpOneLinearIdx = indexerU.linearIndex(i+1,j);
                ssize_t uIpOneJmOneLinearIdx = indexerV.linearIndex(i+1,j-1);
                ssize_t vIpOneLinearIdx = indexerV.linearIndex(i+1,j);

                if(uIpOneLinearIdx != -1
                    && uIpOneJmOneLinearIdx != -1
                    && vIpOneLinearIdx != -1)
                {
                    float lerpedViscosity = viscosityGrid.interpolateAt(fi+0.5f,fj-0.5f);
                    //lerpedViscosity = tempVisc;

                    output.data()[vBaseIndex + idxV]
                                    uIpOneLinearIdx) +=
                        -scaleTwoDx * lerpedViscosity;

                    output.data()[vBaseIndex + idxV]
                                    uIpOneJmOneLinearIdx) +=
                        scaleTwoDx * lerpedViscosity;

                    output.data()[vBaseIndex + idxV]
                                    vBaseIndex + vIpOneLinearIdx) +=
                        -scaleTwoDx * lerpedViscosity;

                    output.data()[vBaseIndex + idxV]
                                    vBaseIndex + idxV) +=
                        scaleTwoDx * lerpedViscosity;
                }
            }
        }
    }

    output.makeCompressed();

    return output; //TEMP!
}

void HeavyViscosityModel::fillRhs(std::vector<double>& rhs, const StaggeredVelocityGrid &velocityGrid, float density)
{
    const LinearIndexable2d& uIndexer = velocityGrid.velocityGridU();
    const LinearIndexable2d& vIndexer = velocityGrid.velocityGridV();
    ssize_t vBaseIndex = uIndexer.linearSize();

    rhs.resize(uIndexer.linearSize() +
               vIndexer.linearSize());

    for (ssize_t i = 0; i < uIndexer.sizeI(); i++)
    {
        for (ssize_t j = 0; j < uIndexer.sizeJ(); j++)
        {
            ssize_t idxU = uIndexer.linearIndex(i,j);

            float u = velocityGrid.getU(i,j);
            rhs[idxU] = density * u;
        }
    }

    for (ssize_t i = 0; i < vIndexer.sizeI(); i++)
    {
        for (ssize_t j = 0; j < vIndexer.sizeJ(); j++)
        {
            ssize_t idxV = vIndexer.linearIndex(i,j);

            float v = velocityGrid.getV(i,j);
            rhs[vBaseIndex + idxV] = density * v;
        }
    }
}

void HeavyViscosityModel::applyResult(StaggeredVelocityGrid &velocityGrid, const std::vector<double> &result)
{
    const LinearIndexable2d& uIndexer = velocityGrid.velocityGridU();
    const LinearIndexable2d& vIndexer = velocityGrid.velocityGridV();
    ssize_t vBaseIndex = velocityGrid.velocityGridU().linearSize();

    for (ssize_t i = 0; i < uIndexer.sizeI(); i++)
    {
        for (ssize_t j = 0; j < uIndexer.sizeJ(); j++)
        {
            ssize_t idxU = uIndexer.linearIndex(i,j);
            velocityGrid.setU(i,j,result[idxU]);
        }
    }

    for (ssize_t i = 0; i < vIndexer.sizeI(); i++)
    {
        for (ssize_t j = 0; j < vIndexer.sizeJ(); j++)
        {
            ssize_t idxV = vIndexer.linearIndex(i,j);
            velocityGrid.setV(i,j,result[vBaseIndex + idxV]);
        }
    }

    if(anyNanInf(velocityGrid.velocityGridU().data()))
    {
        std::cout << "NaN or inf in U velocity after viscosity!\n" << std::flush;
    }

    if(anyNanInf(velocityGrid.velocityGridV().data()))
    {
        std::cout << "NaN or inf in V velocity after viscosity!\n" << std::flush;
    }
}
