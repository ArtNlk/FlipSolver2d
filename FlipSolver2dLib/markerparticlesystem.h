#ifndef MARKERPARTICLESYSTEM_H
#define MARKERPARTICLESYSTEM_H

#include "geometry2d.h"
#include "grid2d.h"
#include "index2d.h"
#include "linearindexable2d.h"
#include "materialgrid.h"
#include "threadpool.h"
#include <array>
#include <variant>
#include <vector>

struct MarkerParticle
{
    Vertex position;
    Vertex velocity;
    float viscosity;
    float temperature;
    float smokeConcentrartion;
    float fuel;
    FluidMaterial material = FluidMaterial::FLUID;
    float testValue = 0.f;
    bool markForDeath = false;
};

template<class T>
inline void swapErase(std::vector<T>& v, size_t index)
{
    std::swap<T>(v[index],v[v.size()-1]);
    v.pop_back();
}

inline void swapErase(std::vector<bool>& v, size_t index)
{
    std::vector<bool>::swap(v[index],v[v.size()-1]);
    v.pop_back();
}

template<class T>
inline T swapPop(std::vector<T>& v, size_t index)
{
    std::swap<T>(v[index],v[v.size()-1]);
    T val = v.back();
    v.pop_back();
    return val;
}

inline bool swapPop(std::vector<bool>& v, size_t index)
{
    std::vector<bool>::swap(v[index],v[v.size()-1]);
    bool val = v.back();
    v.pop_back();
    return val;
}

class RebinSet;

class ParticleBin
{
public:
    ParticleBin(size_t binIdx = std::numeric_limits<size_t>::max());

    using VariantVector = std::variant<std::vector<float>>;

    size_t addMarkerParticle(Vertex position, Vertex velocity = Vertex());

    void eraseMarkerParticle(size_t index);

    void markForDeath(size_t index);

    bool markedForDeath(size_t index);

    Vertex& particlePosition(size_t index);

    Vertex& particleVelocity(size_t index);

    void pruneParticles();

    void rebinParticles(RebinSet &rebinningSet);

    void moveParticleToBin(ParticleBin& other, size_t idx);

    void copyParticleToBin(ParticleBin& other, size_t idx);

    void scheduleRebin(RebinSet &rebinningSet, size_t particleIndex, size_t newBinIndex);

    void clear();

    void setBinIndex(size_t index);

    template<class T>
    std::vector<T>& particleProperties(size_t propertyIndex)
    {
        return std::get<std::vector<T>>(m_properties.at(propertyIndex));
    }

    template<class T>
    size_t addParticleProperty()
    {
        VariantVector v = std::vector<T>(m_particlePositions.size());
        m_properties.push_back(v);
        return m_properties.size() - 1;
    }

    size_t size() const
    {
        return m_particlePositions.size();
    }

    size_t binIdx() const
    {
        return m_binIdx;
    }

    std::vector<Vertex>& positions()
    {
        return m_particlePositions;
    }

    std::vector<Vertex>& velocities()
    {
        return m_velocities;
    }

    std::vector<VariantVector>& properties()
    {
        return m_properties;
    }

protected:

    std::vector<Vertex> m_particlePositions;
    std::vector<Vertex> m_velocities;
    std::vector<bool> m_markedForDeath;
    std::vector<bool> m_markedForRebin;
    std::vector<VariantVector> m_properties;
    size_t m_binIdx;
};

class RebinSet : public ParticleBin
{
public:
    std::vector<size_t>& binIndexes()
    {
        return m_binIndexes;
    }

    size_t particleBinIdx(size_t particleIdx)
    {
        return m_binIndexes.at(particleIdx);
    }

    void clear()
    {
        ParticleBin::clear();

        m_binIndexes.clear();
    }

protected:
    std::vector<size_t> m_binIndexes;
};

enum ParticleAttributeType : size_t
{
    ATTR_FLOAT,
    ATTR_INVALID = std::variant_npos
};

struct RebinRecord
{
    size_t particleIdx;
    size_t newBinIdx;
};

class MarkerParticleSystem
{
public:
    MarkerParticleSystem(int gridSizeI, int gridSizeJ, size_t binSize);

    ParticleBin& binForGridIdx(int linIdx);
    ParticleBin& binForGridIdx(Index2d idx);
    ParticleBin& binForGridIdx(int i, int j);
    ParticleBin& binForGridPosition(Vertex pos);

    ParticleBin& binForBinIdx(int linIdx);

    int gridToBinIdx(Index2d idx);
    int gridToBinIdx(Vertex pos);
    int gridToBinIdx(int linIdx);
    int gridToBinIdx(int i, int j);

    Grid2d<ParticleBin>& bins();

    template<class T>
    size_t addParticleProperty()
    {
        for(auto& bin : m_particleBins.data())
        {
            bin.addParticleProperty<T>();
        }
        m_rebinningSet.addParticleProperty<T>();

        return m_rebinningSet.properties().size() - 1;
    }

    void addMarkerParticle(size_t binIdx, Vertex position, Vertex velocity);

    void markForDeath(size_t binIdx, size_t particleIdx);

    void rebinParticles();

    void scheduleRebin(size_t binIdx, RebinRecord r);

    std::vector<Vertex>& positions(size_t binIdx);

    std::vector<Vertex>& velocities(size_t binIdx);

    Index2d binIdxForIdx(Index2d idx);

    Index2d binIdxForIdx(int i, int j);

    std::array<int,9> binsForGridCell(Index2d idx);

    std::array<int,9> binsForGridCell(int i, int j);

    size_t particleCount() const;

protected:

    void rebinParticlesThread(Range r);

    size_t m_binSize;
    LinearIndexable2d m_gridIndexer;
    Grid2d<ParticleBin> m_particleBins;

    RebinSet m_rebinningSet;
};

#endif // MARKERPARTICLESYSTEM_H
