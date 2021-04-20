#include "PCH.h"
#include "EnergyIntegral.h"
#include "LightingContext.h"
#include "RandomNumberGenerator.h"
#include "CookTorranceBRDF.h"

using namespace DirectX;

static float EstimateEnergy( float cosTheta, float alpha, float etaI, float etaT, uint32_t sampleCount, SRandomNumberGenerator* rng )
{
    XMFLOAT3 wo( sqrtf( 1.0f - cosTheta * cosTheta ), 0.0f, cosTheta );
    float ESum = 0.0f;
    for ( uint32_t iSample = 0; iSample < sampleCount; ++iSample )
    {
        XMFLOAT2 sample = XMFLOAT2( rng->Next(), rng->Next() );

        SLightingContext lightingContext;
        lightingContext.Init( wo );

        XMFLOAT3 wi;
        float EValue;
        float EPdf;
        SampleCookTorranceMicrofacetBRDF( wo, sample, alpha, etaI, etaT, &wi, &EValue, &EPdf, &lightingContext );

        if ( EPdf != 0.0f )
        {
            ESum += EValue * wi.z / EPdf;
        }
    }
    return ESum / sampleCount;
}

void SEnergyIntegral::Execute( uint32_t threadIndex, uint32_t localLaneIndex, uint32_t globalLaneIndex )
{
    /*float cosTheta = localLaneIndex * m_CosThetaInterval;
    float alpha    = ( threadIndex % m_SliceSize ) * m_AlphaInterval;
    float etaT     = m_EtaBegin + ( threadIndex / m_SliceSize ) * m_EtaInterval;
    m_OutputBuffer[ globalLaneIndex ] = EstimateEnergy( cosTheta, alpha, m_EtaI, etaT, m_SampleCount, &m_Rngs[ threadIndex ] );*/
    m_OutputBuffer[ globalLaneIndex ] = globalLaneIndex;
}
