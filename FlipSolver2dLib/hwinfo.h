#ifndef SIMDDISPATCHER_H
#define SIMDDISPATCHER_H

#include <cstdint>
#include <functional>
#include <array>

using SIMDFeatureMask = uint32_t;

enum SIMDFeatures : SIMDFeatureMask {
    SIMD_FEATURE_NONE       =0,
    SIMD_FEATURE_SSE        =(1u << 0),
    SIMD_FEATURE_SSE2       =(1u << 1),
    SIMD_FEATURE_SSE3       =(1u << 2),
    SIMD_FEATURE_SSSE3      =(1u << 3),
    SIMD_FEATURE_SSE41      =(1u << 4),
    SIMD_FEATURE_SSE42      =(1u << 5),
    SIMD_FEATURE_FMA3       =(1u << 6),
    SIMD_FEATURE_AVX        =(1u << 7),
    SIMD_FEATURE_AVX2       =(1u << 8),
    SIMD_FEATURE_AVX512F    =(1u << 9),    //  AVX512 Foundation
    SIMD_FEATURE_AVX512CD   =(1u << 10),   //  AVX512 Conflict Detection
    SIMD_FEATURE_AVX512PF   =(1u << 11),    //  AVX512 Prefetch
    SIMD_FEATURE_AVX512ER   =(1u << 12),    //  AVX512 Exponential + Reciprocal
    SIMD_FEATURE_AVX512VL   =(1u << 13),    //  AVX512 Vector Length Extensions
    SIMD_FEATURE_AVX512BW   =(1u << 14),    //  AVX512 Byte + Word
    SIMD_FEATURE_AVX512DQ   =(1u << 15),    //  AVX512 Doubleword + Quadword
    SIMD_FEATURE_AVX512IFMA =(1u << 16),    //  AVX512 Integer 52-bit Fused Multiply-Add
    SIMD_FEATURE_AVX512VBMI =(1u << 17) //  AVX512 Vector Byte Manipulation Instructions
};

enum SIMDLevel : uint8_t {
    SIMD_LEVEL_NONE = 0,
    SIMD_LEVEL_SSE42,
    SIMD_LEVEL_SSE42_FMA3,
    SIMD_LEVEL_AVX,
    SIMD_LEVEL_AVX2,
    SIMD_LEVEL_AVX512
};

class HwInfo
{
public:
    HwInfo(HwInfo &other) = delete;
    void operator=(const HwInfo &) = delete;

    static HwInfo& i();

    SIMDLevel getSimdLevel();

private:
    HwInfo();
    ~HwInfo() = default;

    SIMDFeatureMask simdFeatures();
    SIMDLevel simdLevel();

    SIMDFeatureMask m_featureMask;
    SIMDLevel m_simdLevel;
};

#endif // SIMDDISPATCHER_H
