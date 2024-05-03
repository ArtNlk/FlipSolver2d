#ifndef PRESSUREDATA_H
#define PRESSUREDATA_H

#include <cstddef>
#include <vector>

#include "customassert.h"
#include "linearindexable2d.h"

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

struct IndexedPressureParameterUnit
{
    IndexedPressureParameterUnit() :
        unitIndex(std::numeric_limits<size_t>::max())
    {
        values[0] = 0.0;
        values[1] = 0.0;
        values[2] = 0.0;
        values[3] = 0.0;
        values[4] = 0.0;
    }

    IndexedPressureParameterUnit(size_t newIdx,
                                 double diag,
                                 double im1,
                                 double ip1,
                                 double jm1,
                                 double jp1) :
        unitIndex(newIdx)
    {
        values[0] = diag;
        values[1] = im1;
        values[2] = ip1;
        values[3] = jm1;
        values[4] = jp1;
    }

    double& diag()
    {
        return values[0];
    }

    double& iNeg()
    {
        return values[1];
    }

    double& iPos()
    {
        return values[2];
    }

    double& jNeg()
    {
        return values[3];
    }

    double& jPos()
    {
        return values[4];
    }

    inline double multiply(const PressureVector& vec,
                           const LinearIndexable2d& indexer) const
    {
        const size_t im1Idx = std::clamp(indexer.linearIdxOfOffset(unitIndex+1,-1,0),size_t(0),vec.trueSize());
        const size_t ip1Idx = std::clamp(indexer.linearIdxOfOffset(unitIndex+1,1,0),size_t(0),vec.trueSize());
        const size_t jm1Idx = std::clamp(indexer.linearIdxOfOffset(unitIndex+1,0,-1),size_t(0),vec.trueSize());
        const size_t jp1Idx = std::clamp(indexer.linearIdxOfOffset(unitIndex+1,0,1),size_t(0),vec.trueSize());
        const size_t centerIdx = unitIndex+1;

        return values[0]*vec.trueAt(centerIdx) +
               values[1]*vec.trueAt(im1Idx) +
               values[2]*vec.trueAt(ip1Idx) +
               values[3]*vec.trueAt(jm1Idx) +
               values[4]*vec.trueAt(jp1Idx);
    }

    inline double multiply(const std::vector<double>& vec,
                           const LinearIndexable2d& indexer) const
    {
        const int im1Idx = indexer.linearIdxOfOffset((int)unitIndex,-1,0);
        const int ip1Idx = indexer.linearIdxOfOffset((int)unitIndex,1,0);
        const int jm1Idx = indexer.linearIdxOfOffset((int)unitIndex,0,-1);
        const int jp1Idx = indexer.linearIdxOfOffset((int)unitIndex,0,1);
        const int centerIdx = (int)unitIndex;

        return values[0]*(centerIdx >= 0 && centerIdx < vec.size() ? vec.at(centerIdx) : 0.0) +
               values[1]*(im1Idx >= 0 && im1Idx < vec.size() ? vec.at(im1Idx) : 0.0) +
               values[2]*(ip1Idx >= 0 && ip1Idx < vec.size() ? vec.at(ip1Idx) : 0.0) +
               values[3]*(jm1Idx >= 0 && jm1Idx < vec.size() ? vec.at(jm1Idx) : 0.0) +
               values[4]*(jp1Idx >= 0 && jp1Idx < vec.size() ? vec.at(jp1Idx) : 0.0);
    }

    size_t unitIndex;
    std::array<double,5> values;
};

class IndexedPressureParameters
{
public:
    IndexedPressureParameters(size_t size, LinearIndexable2d& indexer) :
    m_indexer(indexer)
    {
        m_data.reserve(size);
    }

    void add(IndexedPressureParameterUnit& unit)
    {
        m_data.push_back(unit);
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

        //size_t nextVectorIndex = m_data[0].unitIndex;
        size_t nextDataIndex = 0;
        for(size_t idx = 0; idx < in.size(); idx++)
        {
            if(nextDataIndex < m_data.size() && m_data[nextDataIndex].unitIndex == idx)
            {
                out[idx] = m_data[nextDataIndex].multiply(in,m_indexer);
                nextDataIndex++;
            }
            else
            {
                out[idx] = in[idx];
            }
        }
    }

protected:
    std::vector<IndexedPressureParameterUnit> m_data;
    LinearIndexable2d m_indexer;
};

#endif // PRESSUREDATA_H
