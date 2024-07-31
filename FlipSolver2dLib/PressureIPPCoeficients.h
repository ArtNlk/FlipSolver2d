#ifndef PRESSUREIPPCOEFICIENTS_H
#define PRESSUREIPPCOEFICIENTS_H

#include "linearindexable2d.h"
#include "threadpool.h"

#include <cstddef>
#include <limits>
#include <vector>

struct IndexedIPPCoefficientUnit
{
    IndexedIPPCoefficientUnit() :
        unitIndex(std::numeric_limits<size_t>::max()),
        jNeg(0.0),
        jNegP1(0.0),
        iNeg(0.0),
        iPos(0.0),
        jPosM1(0.0),
        jPos(0.0)
    {

    }

    inline double multiply(const std::vector<double>& vec,
                           const LinearIndexable2d& indexer) const
    {
        const int jNegLinIdx = indexer.linearIdxOfOffset((int)unitIndex,0,-1);
        const int jNegP1LinIdx = indexer.linearIdxOfOffset((int)unitIndex,0,-1)+1;
        const int iNegLinIdx = indexer.linearIdxOfOffset((int)unitIndex,-1,0);
        const int iPosLinIdx = indexer.linearIdxOfOffset((int)unitIndex,1,0);
        const int jPosM1LinIdx = indexer.linearIdxOfOffset((int)unitIndex,0,1)-1;
        const int jPosLinIdx = indexer.linearIdxOfOffset((int)unitIndex,0,1);
        const int centerIdx = (int)unitIndex;

        return  diag   *(centerIdx >= 0 && centerIdx < vec.size() ? vec.at(centerIdx) : 0.0) +
                jNeg   *(jNegLinIdx >= 0 && jNegLinIdx < vec.size() ? vec.at(jNegLinIdx) : 0.0) +
                jNegP1 *(jNegP1LinIdx >= 0 && jNegP1LinIdx < vec.size() ? vec.at(jNegP1LinIdx) : 0.0) +
                iNeg   *(iNegLinIdx >= 0 && iNegLinIdx < vec.size() ? vec.at(iNegLinIdx) : 0.0) +
                iPos   *(iPosLinIdx >= 0 && iPosLinIdx < vec.size() ? vec.at(iPosLinIdx) : 0.0) +
                jPosM1 *(jPosM1LinIdx >= 0 && jPosM1LinIdx < vec.size() ? vec.at(jPosM1LinIdx) : 0.0) +
                jPos   *(jPosLinIdx >= 0 && jPosLinIdx < vec.size() ? vec.at(jPosLinIdx) : 0.0);

    }

    size_t unitIndex;
    double jNeg;
    double jNegP1;
    double iNeg;
    double diag;
    double iPos;
    double jPosM1;
    double jPos;
};

class IndexedIPPCoefficients
{
public:
    IndexedIPPCoefficients(size_t size, LinearIndexable2d& indexer) :
        m_indexer(indexer)
    {
        m_data.reserve(size);
    }

    void add(IndexedIPPCoefficientUnit& unit)
    {
        m_data.push_back(unit);
    }

    void endThreadDataRange()
    {
        if(m_threadDataRanges.empty())
        {
            m_threadDataRanges.push_back(Range(0,m_data.size()));
            return;
        }

        m_threadDataRanges.push_back(Range(m_threadDataRanges.back().end,m_data.size()));
    }

    auto& data()
    {
        return m_data;
    }

    void multiply(const std::vector<double>& in, std::vector<double>& out) const
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

        for(int i = 0; i < ranges.size(); i++)
        {
            if(i >= m_threadDataRanges.size())
            {
                ThreadPool::i()->enqueue(&IndexedIPPCoefficients::multiplyThread,this,
                                         ranges.at(i),Range(0,0),std::cref(in),std::ref(out));
                continue;
            }
            ThreadPool::i()->enqueue(&IndexedIPPCoefficients::multiplyThread,this,
                                     ranges.at(i),m_threadDataRanges.at(i),std::cref(in),std::ref(out));
        }
        ThreadPool::i()->wait();
    }

protected:
    void multiplyThread(Range vecRange, Range dataRange, const std::vector<double>& in, std::vector<double>& out) const
    {
        if(dataRange.size() == 0)
        {
            std::copy(in.begin() + vecRange.start, in.begin() + vecRange.end,out.begin() + vecRange.start);
            return;
        }

        //size_t nextVectorIndex = m_data[0].unitIndex;
        size_t nextDataIndex = dataRange.start;
        for(size_t idx = vecRange.start; idx < vecRange.end; idx++)
        {
            if(nextDataIndex < m_data.size() && nextDataIndex < dataRange.end
                && m_data[nextDataIndex].unitIndex == idx)
            {
                out[idx] = m_data[nextDataIndex].multiply(in, m_indexer);
                nextDataIndex++;
            }
            else
            {
                out[idx] = in[idx];
            }
        }
    }

    std::vector<IndexedIPPCoefficientUnit> m_data;
    std::vector<Range> m_threadDataRanges;
    LinearIndexable2d m_indexer;
};

#endif // PRESSUREIPPCOEFICIENTS_H
