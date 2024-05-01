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

    size_t unitIndex;
    std::array<double,5> values;
};

class IndexedPressureParameters
{
public:
    IndexedPressureParameters(size_t size)
    {
        m_data.reserve(size);
    }

    void add(IndexedPressureParameterUnit& unit)
    {
        m_data.push_back(unit);
    }

protected:
    std::vector<IndexedPressureParameterUnit> m_data;
};

#endif // PRESSUREDATA_H
