#include <chrono>
#include <iostream>

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <ostream>
#include <random>
#include <emmintrin.h>
#include <smmintrin.h>
#include <stdexcept>
#include <thread>

#ifdef FLUID_AVX2
#include <immintrin.h>
#endif

const int gridSize = 4096;
const int innerRepCount = 100000;

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
#ifdef FLUID_AVX2
void matmulAVX2(const std::vector<MatCoeffs>& mat, const std::vector<double>& in, std::vector<double>& out);
#endif

int main()
{
    Eigen::initParallel();
    Eigen::setNbThreads(std::thread::hardware_concurrency());
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
    std::cout << "Eigen running on " << Eigen::nbThreads() << " threads" << std::endl;
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
    mat.makeCompressed();

    auto t1 = high_resolution_clock::now();
    for(int i = 0; i < innerRepCount; i++)
    {
        out = mat*v;
        out.eval();
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

    std::vector<double> outScalar;
    outScalar.resize(gridSize);
    std::vector<double> outSSE;
    outSSE.resize(gridSize);
#ifdef FLUID_AVX2
    std::vector<double> outAVX2;
    outAVX2.resize(gridSize);
#endif

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
    for(int i = 0; i < innerRepCount; i++)
    {
        matmulScalar(mat,v,outScalar);
    }
    auto t2 = high_resolution_clock::now();
    duration<double, std::milli> ms_double = t2 - t1;
    std::cout << "Matmul scalar: "  << ms_double.count() << std::endl;

//    t1 = high_resolution_clock::now();
//    for(int i = 0; i < innerRepCount; i++)
//    {
//        matmulSSE(mat,v,outSSE);
//    }
//    t2 = high_resolution_clock::now();
//    ms_double = t2 - t1;
//    std::cout << "Matmul SSE: "  << ms_double.count() << std::endl;

#ifdef FLUID_AVX2
    t1 = high_resolution_clock::now();
    for(int i = 0; i < innerRepCount; i++)
    {
        matmulAVX2(mat,v,outAVX2);
    }
    t2 = high_resolution_clock::now();
    ms_double = t2 - t1;
    std::cout << "Matmul AVX2: "  << ms_double.count() << std::endl;
#endif

#if 0
    std::cout << "Verifying custom results" << std::endl;

    for(int i = 0; i < outScalar.size(); i++)
    {
        if(outScalar[i] - outSSE[i] > 0.000001)
        {
            std::cout << "Mismatch at "<< i <<"! Scalar vs SSE" << outScalar[i] << " " << outSSE[i] << std::endl;
            throw std::runtime_error("Scalar-SSE mismatch!");
        }
    }
#endif
}

void matmulScalar(const std::vector<MatCoeffs>& mat, const std::vector<double>& in, std::vector<double>& out)
{
#pragma omp parallel for
    for(int rowIdx = 0; rowIdx < mat.size(); rowIdx++)
    {
        const int idxIm1 = std::clamp(rowIdx-1,0,static_cast<int>(in.size())-1);
        const int idxIp1 = std::clamp(rowIdx+1,0,static_cast<int>(in.size())-1);
        const int idxJm1 = std::clamp(rowIdx-gridSize,0,static_cast<int>(in.size())-1);
        const int idxJp1 = std::clamp(rowIdx+gridSize,0,static_cast<int>(in.size())-1);
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
        const int idxIm1 = std::clamp(rowIdx-1,0,static_cast<int>(in.size())-1);
        const int idxIp1 = std::clamp(rowIdx+1,0,static_cast<int>(in.size())-1);
        const int idxJm1 = std::clamp(rowIdx-gridSize,0,static_cast<int>(in.size())-1);
        const int idxJp1 = std::clamp(rowIdx+gridSize,0,static_cast<int>(in.size())-1);
        const double valIdxIm1 = idxIm1 >= 0 ? in[idxIm1] : 0.0;
        const double valIdxIp1 = idxIp1 >= 0 ? in[idxIp1] : 0.0;
        const double valIdxJm1 = idxJm1 < in.size() ? in[idxJm1] : 0.0;
        const double valIdxJp1 = idxJp1 < in.size() ? in[idxJp1] : 0.0;
        const std::array<double,5>& m = mat[rowIdx].vals;
        __m128d _vI, _vJ, _mI, _mJ, _temp, _res;
        _vI = _mm_set_pd(valIdxIm1, valIdxIp1);
        _mI = _mm_loadu_pd(&m[1]);
        _res = _mm_dp_pd(_vI,_mI,0b11001100);

        _vJ = _mm_set_pd(valIdxJm1, valIdxJp1);
        _mJ = _mm_loadu_pd(&m[3]);
        _temp = _mm_dp_pd(_vJ,_mJ,0b11001100);

        out[rowIdx] = m[0]*in[rowIdx];

        _res = _mm_shuffle_pd(_res,_temp,0b01);
        _res = _mm_hadd_pd(_res,_res);

        double result[2] = {0.0,0.0};
        _mm_store_pd(&(result[0]),_res);

        out[rowIdx] += result[0];
//        out[rowIdx] = m[0]*in[rowIdx] +
//                      m[1]*in[idxIm1] +
//                      m[2]*in[idxIp1] +
//                      m[3]*in[idxJm1] +
//                      m[4]*in[idxJp1];
    }
}

#ifdef FLUID_AVX2
void matmulAVX2(const std::vector<MatCoeffs>& mat, const std::vector<double>& in, std::vector<double>& out)
{
#pragma omp parallel for
    for(int rowIdx = 0; rowIdx < mat.size(); rowIdx++)
    {
        const int idxIm1 = std::clamp(rowIdx-1,0,static_cast<int>(in.size()));
        const int idxIp1 = std::clamp(rowIdx+1,0,static_cast<int>(in.size()));
        const int idxJm1 = std::clamp(rowIdx-gridSize,0,static_cast<int>(in.size()));
        const int idxJp1 = std::clamp(rowIdx+gridSize,0,static_cast<int>(in.size()));
        const std::array<double,5>& m = mat[rowIdx].vals;
        __m256d _val, _mat, _temp;
        __m256i _idx;

        _idx = _mm256_set_epi64x(idxJp1, idxJm1, idxIp1, idxIm1);
        _val = _mm256_i64gather_pd(in.data(),_idx,1);
        _mat = _mm256_loadu_pd(m.data());
        _temp = _mm256_mul_pd(_val, _mat);

        __m256d temp = _mm256_hadd_pd( _temp, _temp );
        __m128d hi128 = _mm256_extractf128_pd( temp, 1 );
        __m128d dotproduct = _mm_add_pd( _mm256_castpd256_pd128(temp), hi128 );

        double result[2] = {0.0,0.0};
        _mm_store_pd(&(result[0]),dotproduct);

        out[rowIdx] += result[0];

        //        out[rowIdx] = m[0]*in[rowIdx] +
        //                      m[1]*in[idxIm1] +
        //                      m[2]*in[idxIp1] +
        //                      m[3]*in[idxJm1] +
        //                      m[4]*in[idxJp1];
    }
}
#endif
