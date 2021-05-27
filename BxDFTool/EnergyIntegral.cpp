#include "PCH.h"
#include "EnergyIntegral.h"
#include "RandomNumberGenerator.h"
#include "Utility.h"

using namespace DirectX;
#include "CookTorranceBSDF.h"

template <typename BxDF>
static float EstimateEnergy( float cosTheta, float alpha, float etaI, float etaT, uint32_t sampleCount, SRandomNumberGenerator* rng )
{
    XMFLOAT3 wo( sqrtf( 1.0f - cosTheta * cosTheta ), 0.0f, cosTheta );
    float ESum = 0.0f;
    for ( uint32_t iSample = 0; iSample < sampleCount; ++iSample )
    {
        float selectionSample = BxDF::s_NeedSelectionSample ? rng->Next() : 0.0f;
        XMFLOAT2 bxdfSample = XMFLOAT2( rng->Next(), rng->Next() );

        SLightingContext lightingContext;
        lightingContext.Init( wo );

        XMFLOAT3 wi;
        float value;
        float pdf;
        bool isDeltaBxdf;
        BxDF::Sample( wo, selectionSample, bxdfSample, alpha, etaI, etaT, &wi, &value, &pdf, &isDeltaBxdf, &lightingContext );

        assert( value >= 0.0f && pdf >= 0.0f );

        if ( pdf != 0.0f )
            ESum += value * abs( wi.z ) / pdf;
    }
    return ESum / sampleCount;
}

template <typename BxDF>
void SEnergyIntegral::Execute( uint32_t threadIndex, uint32_t localLaneIndex, uint32_t globalLaneIndex )
{
    float cosTheta = std::max( localLaneIndex * m_CosThetaInterval, 0.0001f );
    float alpha    = ( threadIndex % m_AlphaCount ) * m_AlphaInterval;
    float etaT     = m_EtaBegin + ( threadIndex / m_AlphaCount ) * m_EtaInterval;
    float etaI     = m_EtaI;
    if ( m_Invert ) std::swap( etaI, etaT );
    m_OutputBuffer[ globalLaneIndex ] = EstimateEnergy<BxDF>( cosTheta, alpha, etaI, etaT, m_SampleCount, &m_Rngs[ threadIndex ] );
}

void SEnergyIntegral::Execute_CookTorranceMicrofacetBRDF( uint32_t threadIndex, uint32_t localLaneIndex, uint32_t globalLaneIndex )
{
    Execute<CookTorranceMicrofacetBRDF>( threadIndex, localLaneIndex, globalLaneIndex );
}

void SEnergyIntegral::Execute_CookTorranceMicrofacetBSDF( uint32_t threadIndex, uint32_t localLaneIndex, uint32_t globalLaneIndex )
{
    Execute<CookTorranceMicrofacetBSDF>( threadIndex, localLaneIndex, globalLaneIndex );
}

void SAverageEnergyIntegral::Execute( uint32_t threadIndex, uint32_t localLaneIndex, uint32_t globalLaneIndex )
{
    float cosThetaInterval = 1.0f / ( m_CosThetaCount - 1 );
    float* samples = new float[ m_CosThetaCount ];
    for ( uint32_t i = 0; i < m_CosThetaCount; ++i )
    {
        float cosTheta = std::max( i * cosThetaInterval, 0.0001f );
        samples[ i ] = m_EnergyBuffer[ threadIndex * m_CosThetaCount + i ] * cosTheta;
    }
    float integral = 2.0f * NumericalIntegration::CompositeTrapezoidal( 0.0f, 1.0f, samples, m_CosThetaCount );
    delete[] samples;
    m_OutputBuffer[ threadIndex ] = integral;
}

void SComputeInvCDF::Execute( uint32_t threadIndex, uint32_t localLaneIndex, uint32_t globalLaneIndex )
{
    float cosThetaInterval = 1.0f / ( m_CosThetaCount - 1 );
    float* samples = new float[ m_CosThetaCount ];
    float marginalPDFIntegral = 0.0f;
    float lastMarginalPDFSample = 0.0f;
    for ( uint32_t i = 0; i < m_CosThetaCount; ++i )
    {
        float cosTheta = std::max( i * cosThetaInterval, 0.0001f );
        float E = m_EnergyBuffer[ threadIndex * m_CosThetaCount + i ];

        // Marginal PDF
        float marginalPDFSample = 2.0f * (float)M_PI * std::max( 1.0f - E, 0.0f ) * cosTheta;
        // Calculate the area of the trapezoid
        float area = ( i == 0 ) ? 0.0f : ( lastMarginalPDFSample + marginalPDFSample ) * cosThetaInterval * 0.5f;
        marginalPDFIntegral += area;
        lastMarginalPDFSample = marginalPDFSample;

        float CDF = marginalPDFIntegral;
        samples[ i ] = CDF;
    }

    // Rescale CDF to let CDF(1) = 1
    float CDFMax = marginalPDFIntegral;
    if ( CDFMax != 0.0f )
    {
        float CDFScale = 1.0f / CDFMax;
        for ( uint32_t i = 0; i < m_CosThetaCount; ++i )
        {
            samples[ i ] *= CDFScale;
        }
    }
    samples[ m_CosThetaCount - 1 ] = 1.0f;

    m_PDFScaleOutputBuffer[ threadIndex ] = CDFMax / m_MaxPDFScale;

    // Get inverse CDF
    float* invCDFSamples = m_OutputBuffer + threadIndex * m_CosThetaCount;
    float interval = 1.0f / ( m_CosThetaCount - 1 );
    float lastCDFValue = 0.0f;
    uint32_t iCosTheta = 0;
    for ( uint32_t i = 0; i < m_CosThetaCount; ++i )
    {
        float x = i * interval;
        for ( ; iCosTheta < m_CosThetaCount && samples[ iCosTheta ] < x; ++iCosTheta )
        {
            lastCDFValue = samples[ iCosTheta ];
        }

        float currentCDFValue = samples[ iCosTheta ];
        float t = 1.0f;
        if ( currentCDFValue != lastCDFValue )
        {
            t = ( x - lastCDFValue ) / ( currentCDFValue - lastCDFValue );
        }
        invCDFSamples[ i ] = ( iCosTheta - 1.0f + t ) * interval;
    }
    
    delete[] samples;
}
