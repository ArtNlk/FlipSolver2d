#include "hwinfo.h"


//Thanks to https://github.com/Mysticial/FeatureDetector
#ifdef _WIN32

//  Windows
#include <intrin.h>
#define cpuid(info, x)    __cpuidex(info, x, 0)

#else

//  GCC Intrinsics
#include <cpuid.h>
void cpuid(int info[4], int InfoType){
    __cpuid_count(InfoType, 0, info[0], info[1], info[2], info[3]);
}

#endif

HwInfo::HwInfo()
{
    m_featureMask = simdFeatures();
    m_simdLevel = simdLevel();
}

HwInfo &HwInfo::i()
{
    static HwInfo instance;
    return instance;
}

SIMDLevel HwInfo::getSimdLevel()
{
    return m_simdLevel;
}

SIMDFeatureMask HwInfo::simdFeatures()
{
    SIMDFeatureMask output = SIMD_FEATURE_NONE;
    int info[4];
    cpuid(info, 0);
    int nIds = info[0];

    cpuid(info, 0x80000000);
    unsigned nExIds = info[0];

    //  Detect Features
    if (nIds >= 0x00000001){
        cpuid(info,0x00000001);
        //if((info[3] & ((int)1 << 23)) != 0) output |= SIMD_MMX;
        if((info[3] & ((int)1 << 25)) != 0) output |= SIMD_FEATURE_SSE;
        if((info[3] & ((int)1 << 26)) != 0) output |= SIMD_FEATURE_SSE2;
        if((info[2] & ((int)1 <<  0)) != 0) output |= SIMD_FEATURE_SSE3;

        if((info[2] & ((int)1 <<  9)) != 0) output |= SIMD_FEATURE_SSSE3;
        if((info[2] & ((int)1 << 19)) != 0) output |= SIMD_FEATURE_SSE41;
        if((info[2] & ((int)1 << 20)) != 0) output |= SIMD_FEATURE_SSE42;
        //if((info[2] & ((int)1 << 25)) != 0) output |= SIMD_AES;

        if((info[2] & ((int)1 << 28)) != 0) output |= SIMD_FEATURE_AVX ;
        if((info[2] & ((int)1 << 12)) != 0) output |= SIMD_FEATURE_FMA3;
    }
    if (nIds >= 0x00000007){
        cpuid(info,0x00000007);
        if((info[1] & ((int)1 <<  5)) != 0) output |= SIMD_FEATURE_AVX2;

        if((info[1] & ((int)1 << 16)) != 0) output |= SIMD_FEATURE_AVX512F;
        if((info[1] & ((int)1 << 28)) != 0) output |= SIMD_FEATURE_AVX512CD;
        if((info[1] & ((int)1 << 26)) != 0) output |= SIMD_FEATURE_AVX512PF;
        if((info[1] & ((int)1 << 27)) != 0) output |= SIMD_FEATURE_AVX512ER;
        if((info[1] & ((int)1 << 31)) != 0) output |= SIMD_FEATURE_AVX512VL;
        if((info[1] & ((int)1 << 30)) != 0) output |= SIMD_FEATURE_AVX512BW;
        if((info[1] & ((int)1 << 17)) != 0) output |= SIMD_FEATURE_AVX512DQ;
        if((info[1] & ((int)1 << 21)) != 0) output |= SIMD_FEATURE_AVX512IFMA;
        if((info[2] & ((int)1 <<  1)) != 0) output |= SIMD_FEATURE_AVX512VBMI;
    }
    if (nExIds >= 0x80000001){
        cpuid(info,0x80000001);
        //if((info[3] & ((int)1 << 29)) != 0) output |= SIMD_x64;
        //if((info[2] & ((int)1 <<  5)) != 0) output |= SIMD_ABM;
        //if((info[2] & ((int)1 <<  6)) != 0) output |= SIMD_FEATURE_SSE4a;
        //if((info[2] & ((int)1 << 16)) != 0) output |= SIMD_FEATURE_FMA4;
        //if((info[2] & ((int)1 << 11)) != 0) output |= SIMD_FEATURE_XOP;
    }

    return output;
}

SIMDLevel HwInfo::simdLevel()
{
    SIMDLevel output = SIMD_LEVEL_NONE;

    for(int i = 0; i <= 5; i++)
    {
        if((m_featureMask & (1u << i)) == 0)
        {
            return output;
        }
    }
    output = SIMD_LEVEL_SSE42;

    for(int i = 6; i <= 6; i++)
    {
        if((m_featureMask & (1u << i)) == 0)
        {
            return output;
        }
    }
    output = SIMD_LEVEL_SSE42_FMA3;

    if((m_featureMask & SIMD_FEATURE_AVX) == 0)
    {
        return output;
    }
    output = SIMD_LEVEL_AVX;

    if((m_featureMask & SIMD_FEATURE_AVX2) == 0)
    {
        return output;
    }
    output = SIMD_LEVEL_AVX2;

    for(int i = 9; i <= 17; i++)
    {
        if((m_featureMask & (1u << i)) == 0)
        {
            return output;
        }
    }
    output = SIMD_LEVEL_AVX512;

    return output;
}
