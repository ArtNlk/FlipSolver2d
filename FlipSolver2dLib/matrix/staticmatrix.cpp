#include "staticmatrix.h"
#include "logger.h"

#include "threadpool.h"

#include <stdexcept>
#include <vector>

StaticMatrix::StaticMatrix()
    : m_indexes(0)
    , m_values(0)
    , m_rowStart(0)
    , m_size(0)
{

}

void StaticMatrix::mulThread(Range range, const std::vector<double>& vin, std::vector<double> &vout) const
{
    for(size_t i = range.start; i < range.end; i++)
    {
        vout[i] = 0;
        for(size_t j = m_rowStart[i]; j < m_rowStart[i+1]; j++)
        {
            vout[i] += m_values[j] * vin[m_indexes[j]];
        }
    }
}

void StaticMatrix::multiply(const std::vector<double> &in, std::vector<double> &out) const
{
    //    for(int i = 0; i < size(); i++)
    //    {
    //        //vout[i] = 0;
    //        for(int j = m_rowStart[i]; j < m_rowStart[i+1]; j++)
    //        {
    //            output[i] += m_values[j].second * v[m_values[j].first];
    //        }
    //    }

    std::vector<Range> ranges = ThreadPool::i()->splitRange(out.size());

    for(const Range& range : ranges)
    {
        ThreadPool::i()->enqueue(&StaticMatrix::mulThread,this,range,
                                 std::ref(in),std::ref(out));
    }
    ThreadPool::i()->wait();
}

std::string StaticMatrix::toString()
{
    std::ostringstream output;

    return output.str();
}

size_t StaticMatrix::size() const
{
    return m_size;
}


