#include <chrono>
#include <iostream>

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <ostream>
#include <random>
#include <emmintrin.h>
#include <smmintrin.h>

const int gridSize = 2048;

struct MatCoeffs
{
    MatCoeffs(double diag, double im1, double ip1, double jm1, double jp1)
    {
        vals[0] = diag;
        vals[1] = im1;
        vals[2] = ip1;
        vals[3] = jm1;
        vals[4] = jp1;
    }
    std::array<double,5> vals;
};

struct MatCoeffsIndexed
{
    int idx;
    std::array<double,5> vals;
};

void matmulEigen();
void matmulCustom();

void matmulScalar(const std::vector<MatCoeffs>& mat, const std::vector<double>& in, std::vector<double>& out);
void matmulSSE(const std::vector<MatCoeffs>& mat, const std::vector<double>& in, std::vector<double>& out);
void matmulAVX2(const std::vector<MatCoeffs>& mat, const std::vector<double>& in, std::vector<double>& out);

int main()
{
    for(int i = 0; i < 10; i++)
    {
        std::cout << "Run " << i << std::endl;
        matmulEigen();
        matmulCustom();
    }
    return 0;
}

void matmulEigen()
{
    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;
    using std::chrono::milliseconds;
    std::default_random_engine gen;
    std::uniform_real_distribution<double> dist(0.0,1.0);
    std::vector<Eigen::Triplet<double> > tripletList;
    Eigen::VectorXd v;
    v.resize(gridSize);

    Eigen::VectorXd out;
    out.resize(gridSize);

    for(int i=0;i<gridSize;++i)
    {
        double diag = 0.0;
        if(i-1 >= 0) diag++, tripletList.push_back(Eigen::Triplet(i,i-1,1.0));
        if(i+1 < gridSize) diag++, tripletList.push_back(Eigen::Triplet(i,i+1,1.0));
        if(i-gridSize >= 0) diag++, tripletList.push_back(Eigen::Triplet(i,i-gridSize,1.0));
        if(i+gridSize < gridSize) diag++, tripletList.push_back(Eigen::Triplet(i,i+gridSize,1.0));
        tripletList.push_back(Eigen::Triplet(i,i,diag));
        v.coeffRef(i) = dist(gen);
    }
    Eigen::SparseMatrix<double,Eigen::RowMajor> mat(gridSize,gridSize);
    mat.setFromTriplets(tripletList.begin(), tripletList.end());

    auto t1 = high_resolution_clock::now();
    for(int i = 0; i < 100000; i++)
    {
        out = (mat*v).eval();
    }
    auto t2 = high_resolution_clock::now();
    duration<double, std::milli> ms_double = t2 - t1;
    std::cout << "Eigen: "  << ms_double.count() << std::endl;
}

void matmulCustom()
{
    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;
    using std::chrono::milliseconds;
    std::default_random_engine gen;
    std::uniform_real_distribution<double> dist(0.0,1.0);

    std::vector<MatCoeffs> mat;
    mat.reserve(gridSize);
    std::vector<double> v;
    v.resize(gridSize);

    std::vector<double> out;
    out.resize(gridSize);

    for(int i=0;i<gridSize;++i)
    {
        double diag = 0.0;
        double im1 = 0.0;
        double ip1 = 0.0;
        double jm1 = 0.0;
        double jp1 = 0.0;
        if(i-1 >= 0) diag++, im1 += 1.0;
        if(i+1 < gridSize) diag++, ip1 += 1.0;
        if(i-gridSize >= 0) diag++, jm1 += 1.0;
        if(i+gridSize < gridSize) diag++, jp1 += 1.0;
        mat.emplace_back(diag,im1,ip1,jm1,jp1);
        v.at(i) = dist(gen);
    }

    auto t1 = high_resolution_clock::now();
    for(int i = 0; i < 100000; i++)
    {
        matmulScalar(mat,v,out);
    }
    auto t2 = high_resolution_clock::now();
    duration<double, std::milli> ms_double = t2 - t1;
    std::cout << "Matmul scalar: "  << ms_double.count() << std::endl;

    t1 = high_resolution_clock::now();
    for(int i = 0; i < 100000; i++)
    {
        matmulSSE(mat,v,out);
    }
    t2 = high_resolution_clock::now();
    ms_double = t2 - t1;
    std::cout << "Matmul SSE: "  << ms_double.count() << std::endl;
}

void matmulScalar(const std::vector<MatCoeffs>& mat, const std::vector<double>& in, std::vector<double>& out)
{
#pragma omp parallel for
    for(int rowIdx = 0; rowIdx < mat.size(); rowIdx++)
    {
        const int idxIm1 = rowIdx-1;
        const int idxIp1 = rowIdx+1;
        const int idxJm1 = rowIdx-gridSize;
        const int idxJp1 = rowIdx+gridSize;
        const std::array<double,5>& m = mat[rowIdx].vals;
        out[rowIdx] = m[0]*in[rowIdx] +
                        m[1]*in[idxIm1] +
                        m[2]*in[idxIp1] +
                        m[3]*in[idxJm1] +
                        m[4]*in[idxJp1];
    }
}

void matmulSSE(const std::vector<MatCoeffs>& mat, const std::vector<double>& in, std::vector<double>& out)
{
#pragma omp parallel for
    for(int rowIdx = 0; rowIdx < mat.size(); rowIdx++)
    {
        const int idxIm1 = rowIdx-1;
        const int idxIp1 = rowIdx+1;
        const int idxJm1 = rowIdx-gridSize;
        const int idxJp1 = rowIdx+gridSize;
        const std::array<double,5>& m = mat[rowIdx].vals;
        __m128d _vI, _vJ, _mI, _mJ, _res;
        _vI = _mm_set_pd(in[idxIm1], in[idxIp1]);
        _mI = _mm_loadu_pd(&m[1]);
        _vI = _mm_mul_pd(_vI,_mI);

        _vJ = _mm_set_pd(in[idxJm1], in[idxJp1]);
        _mJ = _mm_loadu_pd(&m[3]);

        out[rowIdx] = m[0]*in[rowIdx];


        _vJ = _mm_mul_pd(_vJ,_mJ);
        _res = _mm_hadd_pd(_vI,_vJ);
        _res = _mm_add_sd(_res,_mm_shuffle_pd(_res, _res, 1));

        double result;
        _mm_store_pd(&result,_res);

        out[rowIdx] += result;
//        out[rowIdx] = m[0]*in[rowIdx] +
//                      m[1]*in[idxIm1] +
//                      m[2]*in[idxIp1] +
//                      m[3]*in[idxJm1] +
//                      m[4]*in[idxJp1];
    }
}
