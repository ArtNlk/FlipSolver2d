#ifndef MATRIXWEIGHTS_H
#define MATRIXWEIGHTS_H

#include "customassert.h"
#include "linearindexable2d.h"
#include "threadpool.h"
#include <vector>

template<class RowUnit>
class MatrixWeights
{
public:
    MatrixWeights(const LinearIndexable2d& indexer) :
    m_indexer(indexer)
    {

    };

    virtual void multiply(const std::vector<double>& in, std::vector<double>& out) const
    {
        ASSERT(in.size() == out.size());

        if(m_data.empty())
        {
            out = in;
            return;
        }

        std::vector<Range> ranges = ThreadPool::i()->splitRange(in.size());

        // multiplyThread(Range(0,in.size()),Range(0,m_data.size()),in,out);

        // return;

        for(size_t i = 0; i < ranges.size(); i++)
        {
            if(i >= m_threadDataRanges.size())
            {
                ThreadPool::i()->enqueue(&MatrixWeights::multiplyThread,this,
                                         ranges.at(i),Range(0,0),std::cref(in),std::ref(out));
                continue;
            }
            ThreadPool::i()->enqueue(&MatrixWeights::multiplyThread,this,
                                     ranges.at(i),m_threadDataRanges.at(i),std::cref(in),std::ref(out));
        }
        ThreadPool::i()->wait();
    }

    virtual void multiplyThread(Range vecRange, Range dataRange, const std::vector<double>& in, std::vector<double>& out) const = 0;

    void endThreadDataRange()
    {
        if (m_threadDataRanges.empty()) {
            m_threadDataRanges.push_back(Range(0, m_data.size()));
            return;
        }

        m_threadDataRanges.push_back(Range(m_threadDataRanges.back().end, m_data.size()));
    }

    std::vector<RowUnit>& data()
    {
        return m_data;
    }

    const std::vector<RowUnit>& data() const
    {
        return m_data;
    }

    void add(RowUnit& unit)
    {
        m_data.push_back(unit);
    }

protected:
    std::vector<Range> m_threadDataRanges;
    std::vector<RowUnit> m_data;
    LinearIndexable2d m_indexer;
};

#endif // MATRIXWEIGHTS_H
