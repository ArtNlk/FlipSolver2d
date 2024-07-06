#include "viscositymodel.h"

ViscosityModel::MatrixType LightViscosityModel::getMatrix(const LinearIndexable2d &indexerU,
                                                         const LinearIndexable2d &indexerV,
                                                         const Grid2d<float> &viscosityGrid,
                                                         const MaterialGrid& materialGrid,
                                                         const float dt,
                                                         const float dx,
                                                         const float density)
{
    const LinearIndexable2d& indexer = viscosityGrid;
    const size_t size = indexerU.linearSize() * 2;

    Eigen::SparseMatrix<double,Eigen::RowMajor> output = Eigen::SparseMatrix<double>();
    output.resize(size,size);
    output.reserve(Eigen::VectorXi::Constant(size,10));

    const double scale = dt;
    const size_t vBaseIndex = indexerU.linearSize();

    for(int i = 0; i < indexer.sizeI()+1; i++)
    {
        for(int j = 0; j < indexer.sizeJ()+1; j++)
        {
            int idxU = indexerU.linearIndex(i,j);
            int idxV = indexerV.linearIndex(i,j);
            float fi = static_cast<float>(i);
            float fj = static_cast<float>(j);

            if(materialGrid.isSolid(i,j))
            {
                if(idxU != -1)
                {
                    output.coeffRef(idxU,idxU) = 1.0;
                }

                if(idxV != -1)
                {
                    output.coeffRef(vBaseIndex + idxV,vBaseIndex + idxV) = 1.0;
                }
                continue;
            }

            if(idxU != -1)
            {
                double diag = 4.0;
                double ip1Neighbor = scale * viscosityGrid.lerpolateAt(fi+1.0,fj);
                double jp1Neighbor = scale * viscosityGrid.lerpolateAt(fi,fj+1.5);
                double im1Neighbor = scale * viscosityGrid.lerpolateAt(fi-1.0,fj);
                double jm1Neighbor = scale * viscosityGrid.lerpolateAt(fi,fj-0.5);

                diag *= viscosityGrid.lerpolateAt(fi-0.5,fj) * scale;
                diag += 1.;

                output.coeffRef(idxU,idxU) = diag;
                if(indexer.inBounds(i+1,j))
                {
                    output.coeffRef(idxU,indexer.linearIndex(i+1,j)) = ip1Neighbor;
                    output.coeffRef(indexer.linearIndex(i+1,j),idxU) = ip1Neighbor;
                }

                if(indexer.inBounds(i,j+1))
                {
                    output.coeffRef(idxU,indexer.linearIndex(i,j+1)) = jp1Neighbor;
                    output.coeffRef(indexer.linearIndex(i,j+1),idxU) = jp1Neighbor;
                }

                if(indexer.inBounds(i-1,j))
                {
                    output.coeffRef(idxU,indexer.linearIndex(i-1,j)) = im1Neighbor;
                    output.coeffRef(indexer.linearIndex(i-1,j),idxU) = im1Neighbor;
                }

                if(indexer.inBounds(i,j-1))
                {
                    output.coeffRef(idxU,indexer.linearIndex(i,j-1)) = jm1Neighbor;
                    output.coeffRef(indexer.linearIndex(i,j-1),idxU) = jm1Neighbor;
                }
            }

            if(idxV != -1)
            {
                double diag = 4.0;
                double ip1Neighbor = scale * viscosityGrid.lerpolateAt(fi+1.5,fj);
                double jp1Neighbor = scale * viscosityGrid.lerpolateAt(fi,fj+1.0);
                double im1Neighbor = scale * viscosityGrid.lerpolateAt(fi-0.5,fj);
                double jm1Neighbor = scale * viscosityGrid.lerpolateAt(fi,fj-1.0);

                diag *= viscosityGrid.lerpolateAt(fi,fj-0.5) * scale;
                diag += 1.;

                output.coeffRef(vBaseIndex+idxV,vBaseIndex+idxV) = diag;
                if(indexer.inBounds(i+1,j))
                {
                    output.coeffRef(vBaseIndex+idxV,indexer.linearIndex(i+1,j)) = ip1Neighbor;
                    output.coeffRef(indexer.linearIndex(i+1,j),vBaseIndex+idxV) = ip1Neighbor;
                }

                if(indexer.inBounds(i,j+1))
                {
                    output.coeffRef(vBaseIndex+idxV,indexer.linearIndex(i,j+1)) = jp1Neighbor;
                    output.coeffRef(indexer.linearIndex(i,j+1),vBaseIndex+idxV) = jp1Neighbor;
                }

                if(indexer.inBounds(i-1,j))
                {
                    output.coeffRef(vBaseIndex+idxV,indexer.linearIndex(i-1,j)) = im1Neighbor;
                    output.coeffRef(indexer.linearIndex(i-1,j),vBaseIndex+idxV) = im1Neighbor;
                }

                if(indexer.inBounds(i,j-1))
                {
                    output.coeffRef(vBaseIndex+idxV,indexer.linearIndex(i,j-1)) = jm1Neighbor;
                    output.coeffRef(indexer.linearIndex(i,j-1),vBaseIndex+idxV) = jm1Neighbor;
                }
            }
        }
    }

    output.makeCompressed();

    return output;
}

ViscosityModel::MatrixType HeavyViscosityModel::getMatrix(const LinearIndexable2d &indexerU,
                                                          const LinearIndexable2d &indexerV,
                                                          const Grid2d<float> &viscosityGrid,
                                                          const MaterialGrid &materialGrid,
                                                          const float dt,
                                                          const float dx,
                                                          const float density)
{
    const LinearIndexable2d& indexer = viscosityGrid;
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
