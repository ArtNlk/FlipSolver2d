#ifndef MARKERPARTICLESYSTEM_H
#define MARKERPARTICLESYSTEM_H

#include "geometry2d.h"
#include "grid2d.h"
#include "index2d.h"
#include "linearindexable2d.h"
#include "materialgrid.h"
#include "threadpool.h"
#include <array>
#include <iterator>
#include <variant>
#include <vector>
#include <latch>

struct MarkerParticle
{
    MarkerParticle(Vertex& pos, Vertex& vel) :
        position(pos),
        velocity(vel),
        testValue(0.f),
        markForDeath(false)
    {

    }

    Vertex position;
    Vertex velocity;
    float testValue = 0.f;
    bool markForDeath = false;
};

class MarkerParticleSystem
{
public:
    using ParticleBin = std::vector<MarkerParticle>;
    using VariantVector = std::variant<std::vector<float>>;
    enum ParticleAttributeType : size_t
    {
        ATTR_FLOAT,
        ATTR_INVALID = std::variant_npos
    };

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

    void rebinParticles();

    void pruneParticles();

    void addMarkerParticle(Vertex position, Vertex velocity = Vertex());

    Grid2d<ParticleBin>& bins();

    Index2d binIdxForIdx(Index2d idx);

    Index2d binIdxForIdx(int i, int j);

    std::array<int,9> binsForGridCell(Index2d idx);

    std::array<int,9> binsForGridCell(int i, int j);

    std::array<int,9> rebinSetForBinIdx(Index2d idx);

protected:

    void rebinParticlesThread(Range r, std::latch& sync);

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

    size_t m_binSize;
    LinearIndexable2d m_gridIndexer;
    Grid2d<ParticleBin> m_particleBins;
    std::vector<MarkerParticle> m_rebinningSet;
    //std::vector<VariantVector> m_properties;
};

#endif // MARKERPARTICLESYSTEM_H
