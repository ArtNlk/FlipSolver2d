#include "viscositymodel.h"
#include "Eigen/src/Core/Matrix.h"

void LightViscosityModel::apply(StaggeredVelocityGrid& velocityGrid,
                              const Grid2d<float>& viscosityGrid,
                              const MaterialGrid& materialGrid,
                              float dt,
                              float dx,
                              float density)
{
    Eigen::VectorXd rhs;
    Eigen::VectorXd result;
    rhs.resize(viscosityGrid.linearSize());

    result.resize(viscosityGrid.linearSize());

    auto viscosityMatrix = getMatrix(velocityGrid,
                                     viscosityGrid,
                                     materialGrid,
                                     dt,
                                     dx,
                                     density);

    m_viscositySolver.setTolerance(1e-4);
    //m_viscositySolver.setMaxIterations(500);
    m_viscositySolver.compute(viscosityMatrix);
    if(m_viscositySolver.info()!=Eigen::Success) {
        std::cout << "Viscosity solver decomposition failed!\n";
        return;
    }

    fillRhs(rhs,velocityGrid.velocityGridU(),viscosityGrid,density);
    result = m_viscositySolver.solve(rhs);
    if(m_viscositySolver.info()!=Eigen::Success) {
        std::cout << "Viscosity solver U solving failed!\n";
        return;
    }
    std::cout << "Viscosity U done with " << m_viscositySolver.iterations() << " iterations\n" << std::endl;

    applyResult(velocityGrid.velocityGridU(), viscosityGrid, result, density);

    fillRhs(rhs,velocityGrid.velocityGridV(),viscosityGrid,density);
    result = m_viscositySolver.solve(rhs);
    if(m_viscositySolver.info()!=Eigen::Success) {
        std::cout << "Viscosity solver V solving failed!\n";
        return;
    }
    std::cout << "Viscosity V done with " << m_viscositySolver.iterations() << " iterations\n";

    applyResult(velocityGrid.velocityGridV(), viscosityGrid, result, density);

    // if(anyNanInf(velocityGrid.velocityGridU().data()))
    // {
    //     std::cout << "NaN or inf in U velocity after viscosity!\n" << std::flush;
    // }

    // if(anyNanInf(velocityGrid.velocityGridV().data()))
    // {
    //     std::cout << "NaN or inf in V velocity after viscosity!\n" << std::flush;
    // }
}

ViscosityModel::MatrixType LightViscosityModel::getMatrix(StaggeredVelocityGrid &velocityGrid,
                                                          const Grid2d<float> &viscosityGrid,
                                                          const MaterialGrid &materialGrid,
                                                          float dt,
                                                          float dx,
                                                          float density)
{
    const LinearIndexable2d& indexer = viscosityGrid;

    const size_t size = indexer.linearSize();

    Eigen::SparseMatrix<double,Eigen::RowMajor> output = Eigen::SparseMatrix<double>();
    output.resize(size,size);
    output.reserve(Eigen::VectorXi::Constant(size,10));

    const double scale = dt;
    //const double scale = 1.0;

    for(int i = 0; i < indexer.sizeI(); i++)
    {
        for(int j = 0; j < indexer.sizeJ(); j++)
        {
            int idx = indexer.linearIndex(i,j);

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

void LightViscosityModel::fillRhs(Eigen::VectorXd& rhs, const Grid2d<float> &velocityGrid, const LinearIndexable2d &indexer, float density)
{
    for (int i = 0; i < indexer.sizeI(); i++)
    {
        for (int j = 0; j < indexer.sizeJ(); j++)
        {
            int idx = indexer.linearIndex(i,j);
            rhs[idx] = density * velocityGrid.at(i,j);
        }
    }
}

void LightViscosityModel::applyResult(Grid2d<float> &velocityGrid, const LinearIndexable2d &indexer, const Eigen::VectorXd &result, float density)
{
    for (int i = 0; i < indexer.sizeI(); i++)
    {
        for (int j = 0; j < indexer.sizeJ(); j++)
        {
            int idx = indexer.linearIndex(i,j);
            velocityGrid.setAt(i,j,result[idx]/density);
        }
    }
}

void HeavyViscosityModel::apply(StaggeredVelocityGrid& velocityGrid,
                              const Grid2d<float>& viscosityGrid,
                              const MaterialGrid& materialGrid,
                              float dt,
                              float dx,
                              float density)
{
    Eigen::VectorXd rhs = getRhs(velocityGrid, density);
    Eigen::VectorXd result;
    result.resize(rhs.size());

    auto viscosityMatrix = getMatrix(velocityGrid,
                                   viscosityGrid,
                                   materialGrid,
                                   dt,
                                   dx,
                                   density);

    m_viscositySolver.setTolerance(1e-4);
    //m_viscositySolver.setMaxIterations(500);
    m_viscositySolver.compute(viscosityMatrix);
    if(m_viscositySolver.info()!=Eigen::Success) {
        std::cout << "Viscosity solver decomposition failed!\n";
        return;
    }

    result = m_viscositySolver.solve(rhs);
    if(m_viscositySolver.info()!=Eigen::Success) {
        std::cout << "Viscosity solver U solving failed!\n";
        return;
    }
    std::cout << "Viscosity done with " << m_viscositySolver.iterations() << " iterations\n" << std::endl;

    applyResult(velocityGrid, result);
}

ViscosityModel::MatrixType HeavyViscosityModel::getMatrix(StaggeredVelocityGrid &velocityGrid,
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

    Eigen::SparseMatrix<double,Eigen::RowMajor> output = Eigen::SparseMatrix<double>();
    output.resize(size,size);
    output.reserve(Eigen::VectorXi::Constant(size,10));

    const float scaleTwoDt = 2*dt / (dx * dx);
    const float scaleTwoDx = dt / (2 * dx * dx);
    const int vBaseIndex = indexerU.linearSize();

    for(int i = 0; i < indexer.sizeI()+1; i++)
    {
        for(int j = 0; j < indexer.sizeJ()+1; j++)
        {
            int idxU = indexerU.linearIndex(i,j);
            int idxV = indexerV.linearIndex(i,j);

            if(idxU != -1)
            {
                float fi = static_cast<float>(i);
                float fj = static_cast<float>(j);

                //U component
                output.coeffRef(idxU,idxU) += density;

                int uImOneLinearIdx = indexerU.linearIndex(i-1,j);

                if(uImOneLinearIdx != -1)
                {
                    output.coeffRef(idxU,
                                    uImOneLinearIdx) += -scaleTwoDt * viscosityGrid.getAt(i-1,j);

                    output.coeffRef(idxU,
                                    idxU) += scaleTwoDt * viscosityGrid.getAt(i-1,j);
                }

                int uIpOneLinearIdx = indexerU.linearIndex(i+1,j);

                if(uIpOneLinearIdx != -1)
                {
                    output.coeffRef(idxU,
                                    uIpOneLinearIdx) += -scaleTwoDt * viscosityGrid.getAt(i,j);

                    output.coeffRef(idxU,
                                    idxU) += scaleTwoDt * viscosityGrid.getAt(i,j);
                }

                int uJmOneLinearIdx = indexerU.linearIndex(i,j-1);
                int vImOneLinearIdx = indexerV.linearIndex(i-1,j);

                if(uJmOneLinearIdx != -1
                    && idxV != -1
                    && vImOneLinearIdx != -1)
                {
                    float lerpedViscosity = viscosityGrid.interpolateAt(fi-0.5f,fj-0.5f);
                    //lerpedViscosity = tempVisc;
                    output.coeffRef(idxU,
                                    uJmOneLinearIdx) += -scaleTwoDx * lerpedViscosity;

                    output.coeffRef(idxU,
                                    vBaseIndex + idxV) += scaleTwoDx * lerpedViscosity;

                    output.coeffRef(idxU,
                                    vBaseIndex + vImOneLinearIdx) += -scaleTwoDx * lerpedViscosity;

                    output.coeffRef(idxU,idxU) += scaleTwoDx * lerpedViscosity;
                }

                int uJpOneLinearIdx = indexerU.linearIndex(i,j+1);
                int vJpOneLinearIdx = indexerV.linearIndex(i,j+1);
                int vImOneJpOneLinearIdx = indexerV.linearIndex(i-1,j+1);

                if(uJpOneLinearIdx != -1
                    && vJpOneLinearIdx != -1
                    && vImOneJpOneLinearIdx != -1)
                {
                    float lerpedViscosity = viscosityGrid.interpolateAt(fi-0.5f,fj+0.5f);
                    //lerpedViscosity = tempVisc;
                    output.coeffRef(idxU,
                                    uJpOneLinearIdx) += -scaleTwoDx * lerpedViscosity;

                    output.coeffRef(idxU,
                                    vBaseIndex + vJpOneLinearIdx) += -scaleTwoDx * lerpedViscosity;

                    output.coeffRef(idxU,
                                    vBaseIndex + vImOneJpOneLinearIdx) += scaleTwoDx * lerpedViscosity;

                    output.coeffRef(idxU,idxU) += scaleTwoDx * lerpedViscosity;
                }
            }

            //V component
            if(idxV != -1)
            {
                float fi = static_cast<float>(i);
                float fj = static_cast<float>(j);
                int vBaseIndex = indexerU.linearSize();
                output.coeffRef(vBaseIndex + idxV,vBaseIndex + idxV) += density;

                int vJmOneLinearIdx = indexerV.linearIndex(i,j-1);

                if(vJmOneLinearIdx != -1)
                {
                    output.coeffRef(vBaseIndex + idxV,
                                    vBaseIndex + vJmOneLinearIdx) += -scaleTwoDt * viscosityGrid.getAt(i,j-1);

                    output.coeffRef(vBaseIndex + idxV,
                                    vBaseIndex + idxV) += scaleTwoDt * viscosityGrid.getAt(i,j-1);
                }

                int vJpOneLinearIdx = indexerV.linearIndex(i,j+1);

                if(vJpOneLinearIdx != -1)
                {
                    output.coeffRef(vBaseIndex + idxV,
                                    vBaseIndex + vJpOneLinearIdx)
                        += -scaleTwoDt * viscosityGrid.getAt(i,j);

                    output.coeffRef(vBaseIndex + idxV,
                                    vBaseIndex + idxV)
                        += scaleTwoDt * viscosityGrid.getAt(i,j);
                }

                int uJmOneLinearIdx = indexerU.linearIndex(i,j-1);
                int vImOneLinearIdx = indexerV.linearIndex(i-1,j);

                if(idxU != -1
                    && uJmOneLinearIdx != -1
                    && vImOneLinearIdx != -1)
                {
                    float lerpedViscosity = viscosityGrid.interpolateAt(fi-0.5f,fj-0.5f);
                    //lerpedViscosity = tempVisc;

                    output.coeffRef(vBaseIndex + idxV,
                                    idxU) +=
                        scaleTwoDx * lerpedViscosity;

                    output.coeffRef(vBaseIndex + idxV,
                                    uJmOneLinearIdx) +=
                        -scaleTwoDx * lerpedViscosity;

                    output.coeffRef(vBaseIndex + idxV,
                                    vBaseIndex + vImOneLinearIdx) +=
                        -scaleTwoDx * lerpedViscosity;

                    output.coeffRef(vBaseIndex + idxV,
                                    vBaseIndex + idxV) +=
                        scaleTwoDx * lerpedViscosity;
                }

                int uIpOneLinearIdx = indexerU.linearIndex(i+1,j);
                int uIpOneJmOneLinearIdx = indexerV.linearIndex(i+1,j-1);
                int vIpOneLinearIdx = indexerV.linearIndex(i+1,j);

                if(uIpOneLinearIdx != -1
                    && uIpOneJmOneLinearIdx != -1
                    && vIpOneLinearIdx != -1)
                {
                    float lerpedViscosity = viscosityGrid.interpolateAt(fi+0.5f,fj-0.5f);
                    //lerpedViscosity = tempVisc;

                    output.coeffRef(vBaseIndex + idxV,
                                    uIpOneLinearIdx) +=
                        -scaleTwoDx * lerpedViscosity;

                    output.coeffRef(vBaseIndex + idxV,
                                    uIpOneJmOneLinearIdx) +=
                        scaleTwoDx * lerpedViscosity;

                    output.coeffRef(vBaseIndex + idxV,
                                    vBaseIndex + vIpOneLinearIdx) +=
                        -scaleTwoDx * lerpedViscosity;

                    output.coeffRef(vBaseIndex + idxV,
                                    vBaseIndex + idxV) +=
                        scaleTwoDx * lerpedViscosity;
                }
            }
        }
    }

    output.makeCompressed();

    return output;
}

Eigen::VectorXd HeavyViscosityModel::getRhs(const StaggeredVelocityGrid &velocityGrid, float density)
{
    const LinearIndexable2d& uIndexer = velocityGrid.velocityGridU();
    const LinearIndexable2d& vIndexer = velocityGrid.velocityGridV();
    int vBaseIndex = uIndexer.linearSize();

    Eigen::VectorXd rhs;
    rhs.resize(uIndexer.linearSize() +
               vIndexer.linearSize());

    for (int i = 0; i < uIndexer.sizeI(); i++)
    {
        for (int j = 0; j < uIndexer.sizeJ(); j++)
        {
            int idxU = uIndexer.linearIndex(i,j);

            float u = velocityGrid.getU(i,j);
            rhs[idxU] = density * u;
        }
    }

    for (int i = 0; i < vIndexer.sizeI(); i++)
    {
        for (int j = 0; j < vIndexer.sizeJ(); j++)
        {
            int idxV = vIndexer.linearIndex(i,j);

            float v = velocityGrid.getV(i,j);
            rhs[vBaseIndex + idxV] = density * v;
        }
    }

    return rhs;
}

void HeavyViscosityModel::applyResult(StaggeredVelocityGrid &velocityGrid,
                                      const Eigen::VectorXd &result)
{
    const LinearIndexable2d& uIndexer = velocityGrid.velocityGridU();
    const LinearIndexable2d& vIndexer = velocityGrid.velocityGridV();
    int vBaseIndex = velocityGrid.velocityGridU().linearSize();

    for (int i = 0; i < uIndexer.sizeI(); i++)
    {
        for (int j = 0; j < uIndexer.sizeJ(); j++)
        {
            int idxU = uIndexer.linearIndex(i,j);
            velocityGrid.setU(i,j,result[idxU]);
        }
    }

    for (int i = 0; i < vIndexer.sizeI(); i++)
    {
        for (int j = 0; j < vIndexer.sizeJ(); j++)
        {
            int idxV = vIndexer.linearIndex(i,j);
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
