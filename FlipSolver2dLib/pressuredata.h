#ifndef PRESSUREDATA_H
#define PRESSUREDATA_H

#include <cstddef>
#include <cstdint>
#include <vector>

#include "customassert.h"
#include "linearindexable2d.h"
#include "threadpool.h"

class PressureVector
{
public:
    PressureVector(size_t size = 0, double init = 0.0) :
        m_data(size+2, init)
    {

    }

    size_t size() const
    {
        return m_data.size() - 2;
    }

    size_t trueSize() const
    {
        return m_data.size();
    }

    double& operator[](size_t index)
    {
        return m_data[index+1];
    }

    double& at(size_t index)
    {
        ASSERT(index < m_data.size() - 1);
        return m_data.at(index+1);
    }

    double& trueAt(size_t index)
    {
        return m_data.at(index);
    }

    double operator[](size_t index) const
    {
        return m_data[index+1];
    }

    double at(size_t index) const
    {
        ASSERT(index < m_data.size() - 1);
        return m_data.at(index+1);
    }

    double trueAt(size_t index) const
    {
        return m_data.at(index);
    }

    std::vector<double>::const_iterator cbegin() const
    {
        return m_data.cbegin()+1;
    }

    std::vector<double>::const_iterator cend() const
    {
        return m_data.cend()-1;
    }

    std::vector<double>::iterator begin()
    {
        return m_data.begin()+1;
    }

    std::vector<double>::iterator end()
    {
        return m_data.end()-1;
    }

protected:
    std::vector<double> m_data;
};

enum NeighborMask : uint8_t
{
    I_NEG_NEIGHBOR_BIT = 1<<0,
    I_POS_NEIGHBOR_BIT = 1<<1,
    J_NEG_NEIGHBOR_BIT = 1<<2,
    J_POS_NEIGHBOR_BIT = 1<<3
    // K_NEG_NEIGHBOR_BIT = 1<<4
    // K_POS_NEIGHBOR_BIT = 1<<5
};

struct IndexedPressureParameterUnit
{
    IndexedPressureParameterUnit() :
        unitIndex(std::numeric_limits<size_t>::max()),
        fluidNeighborMask(0),
        nonsolidNeighborCount(0)
    {

    }

    void setNeighbor(NeighborMask bit, bool flag)
    {
        if (flag) {
            fluidNeighborMask |= bit;
        } else {
            fluidNeighborMask &= ~bit;
        }
    }

    // inline double multiply(const PressureVector& vec,
    //                        const LinearIndexable2d& indexer) const
    // {
    //     const size_t im1Idx = std::clamp(indexer.linearIdxOfOffset(unitIndex+1,-1,0),size_t(0),vec.trueSize());
    //     const size_t ip1Idx = std::clamp(indexer.linearIdxOfOffset(unitIndex+1,1,0),size_t(0),vec.trueSize());
    //     const size_t jm1Idx = std::clamp(indexer.linearIdxOfOffset(unitIndex+1,0,-1),size_t(0),vec.trueSize());
    //     const size_t jp1Idx = std::clamp(indexer.linearIdxOfOffset(unitIndex+1,0,1),size_t(0),vec.trueSize());
    //     const size_t centerIdx = unitIndex+1;

    //     return values[0]*vec.trueAt(centerIdx) +
    //            values[1]*vec.trueAt(im1Idx) +
    //            values[2]*vec.trueAt(ip1Idx) +
    //            values[3]*vec.trueAt(jm1Idx) +
    //            values[4]*vec.trueAt(jp1Idx);
    // }

    inline double multiply(const std::vector<double>& vec,
                           const LinearIndexable2d& indexer, const double scale) const
    {
        const ssize_t im1Idx = indexer.linearIdxOfOffset(unitIndex,-1,0);
        const ssize_t ip1Idx = indexer.linearIdxOfOffset(unitIndex,1,0);
        const ssize_t jm1Idx = indexer.linearIdxOfOffset(unitIndex,0,-1);
        const ssize_t jp1Idx = indexer.linearIdxOfOffset(unitIndex,0,1);
        const ssize_t centerIdx = unitIndex;

        return (scale * nonsolidNeighborCount)*(centerIdx >= 0 && centerIdx < vec.size() ? vec.at(centerIdx) : 0.0) +
               (-scale*(fluidNeighborMask&I_NEG_NEIGHBOR_BIT))*(im1Idx >= 0 && im1Idx < vec.size() ? vec.at(im1Idx) : 0.0) +
               (-scale*((fluidNeighborMask&I_POS_NEIGHBOR_BIT)>>1))*(ip1Idx >= 0 && ip1Idx < vec.size() ? vec.at(ip1Idx) : 0.0) +
               (-scale*((fluidNeighborMask&J_NEG_NEIGHBOR_BIT)>>2))*(jm1Idx >= 0 && jm1Idx < vec.size() ? vec.at(jm1Idx) : 0.0) +
               (-scale*((fluidNeighborMask&J_POS_NEIGHBOR_BIT)>>3))*(jp1Idx >= 0 && jp1Idx < vec.size() ? vec.at(jp1Idx) : 0.0);
    }

    size_t unitIndex;
    uint8_t fluidNeighborMask;
    uint8_t nonsolidNeighborCount;
};

class IndexedPressureParameters
{
public:
    IndexedPressureParameters(size_t size, LinearIndexable2d& indexer, double scale) :
    m_indexer(indexer),
    m_scale(scale)
    {
        m_data.reserve(size);
    }

    void add(IndexedPressureParameterUnit& unit)
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

    const auto& data() const
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

        for(size_t i = 0; i < ranges.size(); i++)
        {
            if(i >= m_threadDataRanges.size())
            {
                ThreadPool::i()->enqueue(&IndexedPressureParameters::multiplyThread,this,
                                         ranges.at(i),Range(0,0),std::cref(in),std::ref(out));
                continue;
            }
            ThreadPool::i()->enqueue(&IndexedPressureParameters::multiplyThread,this,
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
                out[idx] = m_data[nextDataIndex].multiply(in, m_indexer, m_scale);
                nextDataIndex++;
            }
            else
            {
                out[idx] = in[idx];
            }
        }
    }

    std::vector<IndexedPressureParameterUnit> m_data;
    std::vector<Range> m_threadDataRanges;
    LinearIndexable2d m_indexer;
    double m_scale;
};

#endif // PRESSUREDATA_H
