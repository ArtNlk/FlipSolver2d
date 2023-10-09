#include "vmath.h"
#ifdef FLUID_SSE
#include "vops_sse42.h"
#elif defined FLUID_AVX2
#include "vops_avx2.h"
#endif

void VOps::dotThread(Range range, std::vector<double> &v1, std::vector<double> &v2, double &output)
{
    ASSERT(v1.size() == v2.size());
    double result = 0.0;
    for(int i = range.start; i < range.end; i++)
    {
        result += v1[i] * v2[i];
    }
    output = result;
}

VOps &VOps::i()
{
    static VOps instance;
    return instance;
}

double VOps::dot(std::vector<double> &v1, std::vector<double> &v2) {
    ASSERT(v1.size() == v2.size());

    double result = 0;
    std::vector<Range> ranges = ThreadPool::i()->splitRange(v1.size());
    std::vector<double> outputs(ranges.size());

    int i = 0;
    for (Range &range : ranges) {
        ThreadPool::i()->enqueue(m_dotThread, range, std::ref(v1), std::ref(v2),
                                 std::ref(outputs[i]));
        i++;
    }
    ThreadPool::i()->wait();

    for (double &v : outputs) {
        result += v;
    }

    return result;
}

bool VOps::isZero(const std::vector<double> &v1, const double eps)
{
    for(int i = 0; i < v1.size(); i++)
    {
        if(std::abs(v1[i]) > eps)
        {
            return false;
        }
    }

    return true;
}

void VOps::addMulThread(Range range, std::vector<double> &output, const std::vector<double> &vec1, const std::vector<double> &vec2, double value)
{
    for(int i = range.start; i < range.end; i++)
    {
        output[i] = vec1[i] + vec2[i]*value;
    }
}

void VOps::addMul(std::vector<double> &output, const std::vector<double> &vec1, const std::vector<double> &vec2, double value)
{
    ASSERT(output.size() == vec1.size() && output.size() == vec2.size());
    std::vector<Range> ranges = ThreadPool::i()->splitRange(vec1.size());

    for(Range& range : ranges)
    {
        ThreadPool::i()->enqueue(m_addMulThread,range,std::ref(output),std::ref(vec1),std::ref(vec2),value);
    }
    ThreadPool::i()->wait();
}

void VOps::subMulThread(Range range, std::vector<double> &output, const std::vector<double> &vec1, const std::vector<double> &vec2, double value)
{
    for(int i = range.start; i < range.end; i++)
    {
        output[i] = vec1[i] - vec2[i]*value;
    }
}

void VOps::subMul(std::vector<double> &output, const std::vector<double> &vec1, const std::vector<double> &vec2, double value)
{
    ASSERT(output.size() == vec1.size() && output.size() == vec2.size());
    std::vector<Range> ranges = ThreadPool::i()->splitRange(vec1.size());

    for(Range& range : ranges)
    {
        ThreadPool::i()->enqueue(m_subMulThread,range,std::ref(output),std::ref(vec1),std::ref(vec2),value);
    }
    ThreadPool::i()->wait();
}

void VOps::maxAbsThread(Range range, std::vector<double> &vec, double &max)
{
    double maxAbs = std::numeric_limits<double>::min();
    for(int i = range.start; i < range.end; i++)
    {
        if(std::abs(vec[i]) > max)
        {
            maxAbs = std::abs(vec[i]);
        }
    }
    max = maxAbs;
}

double VOps::maxAbs(std::vector<double> &vec)
{
    std::vector<Range> ranges = ThreadPool::i()->splitRange(vec.size());
    std::vector<double> maxAbs(ranges.size());
    double max = 0;

    int i = 0;
    for(Range& range : ranges)
    {
        ThreadPool::i()->enqueue(m_maxAbsThread,range,std::ref(vec),std::ref(maxAbs[i]));
        i++;
    }
    ThreadPool::i()->wait();

    for(double& v : maxAbs)
    {
        if(v > max)
        {
            max = std::abs(v);
        }
    }

    return max;
}

VOps::VOps() :
    m_dotThread(&VOps::dotThread),
    m_addMulThread(&VOps::addMulThread),
    m_subMulThread(&VOps::subMulThread),
    m_maxAbsThread(&VOps::maxAbsThread)
{
#ifdef FLUID_SSE
        m_addMulThread = &VOps_sse42::addMulThread;
#elif defined FLUID_AVX
        m_addMulThread = &VOps_avx::addMulThread;
#endif
}
