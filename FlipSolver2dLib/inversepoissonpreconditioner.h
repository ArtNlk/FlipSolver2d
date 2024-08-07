#ifndef INVERSEPOISSONPRECONDITIONER_H
#define INVERSEPOISSONPRECONDITIONER_H

#include <Eigen/Core>
#include <Eigen/Sparse>

template <typename Scalar, int UpLo_ = Eigen::Lower, typename OrderingType_ = Eigen::AMDOrdering<int> >
class InversePoissonPreconditioner : public Eigen::SparseSolverBase<InversePoissonPreconditioner<Scalar,UpLo_,OrderingType_> >
{
protected:
    typedef Eigen::SparseSolverBase<InversePoissonPreconditioner<Scalar,UpLo_,OrderingType_> > Base;
    using Base::m_isInitialized;
public:
    typedef typename Eigen::NumTraits<Scalar>::Real RealScalar;
    typedef OrderingType_ OrderingType;
    typedef typename OrderingType::PermutationType PermutationType;
    typedef typename PermutationType::StorageIndex StorageIndex;
    typedef Eigen::SparseMatrix<Scalar,Eigen::RowMajor,StorageIndex> FactorType;
    enum { UpLo = UpLo_ };
    enum {
        ColsAtCompileTime = Eigen::Dynamic,
        MaxColsAtCompileTime = Eigen::Dynamic
    };
public:

    InversePoissonPreconditioner() : m_analysisIsOk(false),m_factorizationIsOk(false) {}

    template<typename MatrixType>
    InversePoissonPreconditioner(const MatrixType& matrix) : m_analysisIsOk(false),m_factorizationIsOk(false)
    {
        compute(matrix);
    }

    EIGEN_CONSTEXPR Eigen::Index rows() const EIGEN_NOEXCEPT { return m_mat.rows(); }

    EIGEN_CONSTEXPR Eigen::Index cols() const EIGEN_NOEXCEPT { return m_mat.cols(); }


    Eigen::ComputationInfo info() const
    {
        eigen_assert(m_isInitialized && "InversePoissonPreconditioner is not initialized.");
        return m_info;
    }

    void factorize(const Eigen::SparseMatrix<double,Eigen::RowMajor>& mat)
    {
        auto identity = mat;
        identity.setIdentity();
        Eigen::VectorXd invdiag;
        invdiag.resize(mat.cols());

#pragma omp parallel for
        for(ssize_t j=0; j<mat.outerSize(); ++j)
        {
            Eigen::SparseMatrix<double,Eigen::RowMajor>::InnerIterator it(mat,j);
            while(it && it.index()!=j) ++it;
            if(it && it.index()==j && it.value()!=Scalar(0))
                invdiag(j) = Scalar(1)/it.value();
            else
                invdiag(j) = Scalar(1);
        }
        m_mat = mat.template triangularView<Eigen::StrictlyLower>();

#pragma omp parallel for
        for(ssize_t j=0; j<m_mat.outerSize(); ++j)
        {
            Eigen::SparseMatrix<double,Eigen::RowMajor>::InnerIterator it(m_mat,j);
            for(;it && it.index() < j; ++it)
            {
                it.valueRef() *= invdiag(it.col());
            }
        }
        m_mat = identity - m_mat;
        m_mat = m_mat * m_mat.transpose();

        m_info = Eigen::Success;
        m_factorizationIsOk = true;
        return;
    }

    template<typename MatrixType>
    void compute(const MatrixType& mat)
    {
        factorize(mat);
        m_isInitialized = true;
    }

    // internal
    template<typename Rhs, typename Dest>
    void _solve_impl(const Rhs& b, Dest& x) const
    {
        eigen_assert(m_factorizationIsOk && "factorize() should be called first");
        x = m_mat * b;
    }

protected:
    FactorType m_mat;
    bool m_analysisIsOk;
    bool m_factorizationIsOk;
    Eigen::ComputationInfo m_info;
    PermutationType m_perm;
};

#endif // INVERSEPOISSONPRECONDITIONER_H
