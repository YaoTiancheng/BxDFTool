// BRDFIntegral.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "PCH.h"

#include "RandomNumberGenerator.h"
//#include "CookTorranceBRDF.h"
#include "EnergyIntegral.h"
#include "Utility.h"
#include "MultiProcessing/MPDispatcher.h"

using namespace DirectX;

const uint32_t s_RoughnessCount = 32;
const uint32_t s_CosThetaCount = 32;
const uint32_t s_IORCount = 1;
const float    s_RoughnessStep = 1.0f / ( s_RoughnessCount - 1 );

//void CalculateEnergyForAllCosTheta( float alpha, float etaI, float etaT, SRandomNumberGenerator* rng, uint32_t sampleCount, float* EBuffer )
//{
//    const float s_CosThetaStep = 1.0f / ( s_CosThetaCount - 1 );
//    for ( uint32_t iCosTheta = 0; iCosTheta < s_CosThetaCount; ++iCosTheta )
//    {
//        float cosTheta = std::max( iCosTheta * s_CosThetaStep, 0.001f );
//        float EEstimate = EstimateEnergy( cosTheta, alpha, etaI, etaT, sampleCount, rng );
//        EBuffer[ iCosTheta ] = EEstimate;
//    }
//}
//
//float EstimateAverageEnergy( float alpha, float etaI, float etaT, SRandomNumberGenerator* rng, uint32_t cosThetaSampleCount, uint32_t wiSampleCount )
//{
//    float EAvgSum = 0.0f;
//    for ( uint32_t iSample = 0; iSample < cosThetaSampleCount; ++iSample )
//    {
//        float cosTheta = rng->Next();
//        EAvgSum += EstimateEnergy( cosTheta, alpha, etaI, etaT, wiSampleCount, rng ) * cosTheta;
//    }
//    return 2.0f * EAvgSum / cosThetaSampleCount;
//}
//
//float EstimateAverageEnergyTrapezoidal( float alpha, float etaI, float etaT, SRandomNumberGenerator* rng, uint32_t cosThetaSampleCount, uint32_t wiSampleCount )
//{
//    float* samples = new float[ cosThetaSampleCount ];
//    float cosThetaInterval = 1.0f / ( cosThetaSampleCount - 1 );
//    for ( uint32_t iSample = 0; iSample < cosThetaSampleCount; ++iSample )
//    {
//        float cosTheta = std::max( iSample * cosThetaInterval, 0.001f );
//        samples[ iSample ] = EstimateEnergy( cosTheta, alpha, etaI, etaT, wiSampleCount, rng ) * cosTheta;
//    }
//    float estimate = 2.0f * NumericalIntegration::CompositeTrapezoidal( 0.0f, 1.0f, samples, cosThetaSampleCount );
//    delete[] samples;
//    return estimate;
//}
//
//float EstimateAverageFresnel( float alpha, float etaI, float etaT, SRandomNumberGenerator* rng )
//{
//    const uint32_t s_SampleCount = 18000;
//    float FAvgSum = 0.0f;
//    for ( uint32_t iSample = 0; iSample < s_SampleCount; ++iSample )
//    {
//        float cosTheta = rng->Next();
//
//        XMFLOAT3 m;
//        SampleGGXMicrofacetDistribution( XMFLOAT2( rng->Next(), rng->Next() ), alpha, &m );
//
//        XMFLOAT3 wo( sqrtf( 1.0f - cosTheta * cosTheta ), 0.0f, cosTheta );
//
//        float WOdotM;
//        DirectX::XMStoreFloat( &WOdotM, DirectX::XMVector3Dot( DirectX::XMLoadFloat3( &wo ), DirectX::XMLoadFloat3( &m ) ) );
//
//        if ( WOdotM >= 0.0f )
//        {
//            FAvgSum += EvaluateDielectricFresnel( WOdotM, etaI, etaT ) * cosTheta;
//        }
//    }
//    return 2.0f * FAvgSum / s_SampleCount;
//}

void PrintEnergyBuffer( float* buffer, uint32_t column, uint32_t row, uint32_t slice )
{
    float* pBuffer = buffer;
    for ( uint32_t iSlice = 0; iSlice < slice; ++iSlice )
    {
        for ( uint32_t iRow = 0; iRow < row; ++iRow )
        {
            for ( uint32_t iColumn = 0; iColumn < column; ++iColumn )
            {
                float value = *pBuffer;
                //value = std::roundf( std::min( 1.0f, value ) * 65535.0f );
                printf( "%9.7ff, ", *pBuffer );
                //printf( "%5d, ", (uint16_t)value );
                pBuffer++;
            }
            printf( "\n" );
        }
        printf( "\n" );
    }
}

//struct SAverageEnergyEstimationJob
//{
//    uint32_t m_RoughnessBegin;
//    uint32_t m_RoughnessEnd;
//    uint32_t m_CosThetaSampleCount;
//    uint32_t m_WiSampleCount;
//    SRandomNumberGenerator m_Rng;
//    float* m_ResultBuffer;
//
//    void Execute()
//    {
//        for ( uint32_t iRoughness = m_RoughnessBegin; iRoughness < m_RoughnessEnd; ++iRoughness )
//        {
//            float roughness = s_RoughnessStep * iRoughness;
//            m_ResultBuffer[ iRoughness ] = EstimateAverageEnergyTrapezoidal( roughness, 0.0f, 1.0f, &m_Rng, m_CosThetaSampleCount, m_WiSampleCount );
//        }
//    }
//};

int main()
{
    uint32_t estimateTarget = 0;

    if ( estimateTarget == 0 )
    {
        bool calculateE = false;

        float EBuffer[ s_RoughnessCount * s_CosThetaCount * s_IORCount ];
        memset( EBuffer, 0, sizeof( float ) * s_RoughnessCount * s_CosThetaCount * s_IORCount );
        SRandomNumberGenerator rngs[ s_RoughnessCount * s_IORCount ];

        SEnergyIntegral integral;
        integral.m_AlphaInterval = s_RoughnessStep;
        integral.m_CosThetaInterval = 1.0f / ( s_CosThetaCount - 1 );
        integral.m_EtaBegin = 1.0f;
        integral.m_EtaInterval = ( 3.0f - integral.m_EtaBegin ) / ( s_IORCount - 1 );
        integral.m_EtaI = calculateE ? 0.0f : 1.0f;
        integral.m_SliceSize = s_CosThetaCount * s_RoughnessCount;
        integral.m_OutputBuffer = EBuffer;
        integral.m_Rngs = rngs;
        integral.m_SampleCount = 12000;

        using std::placeholders::_1;
        using std::placeholders::_2;
        using std::placeholders::_3;
        MP::LaneFunctionType laneFunction = std::bind( &SEnergyIntegral::Execute, &integral, _1, _2, _3 );
        MP::Dispatch( s_CosThetaCount, s_RoughnessCount * s_IORCount, laneFunction );

        PrintEnergyBuffer( EBuffer, s_CosThetaCount, s_RoughnessCount, s_IORCount );
    }
    else if ( estimateTarget == 1 )
    {
        /*float EAvgBuffer[ s_RoughnessCount ];
        SAverageEnergyEstimationJob jobs[ 4 ];
        for ( auto& job : jobs )
        {
            job.m_ResultBuffer = EAvgBuffer;
            job.m_WiSampleCount = 960000;
            job.m_CosThetaSampleCount = 32;
        }
        jobs[ 0 ].m_RoughnessBegin = 0;
        jobs[ 0 ].m_RoughnessEnd = 8;
        jobs[ 1 ].m_RoughnessBegin = 8;
        jobs[ 1 ].m_RoughnessEnd = 16;
        jobs[ 2 ].m_RoughnessBegin = 16;
        jobs[ 2 ].m_RoughnessEnd = 24;
        jobs[ 3 ].m_RoughnessBegin = 24;
        jobs[ 3 ].m_RoughnessEnd = 32;
        std::thread* threads[ 4 ];
        for ( uint32_t i = 0; i < 4; ++i )
        {
            threads[ i ] = new std::thread( &SAverageEnergyEstimationJob::Execute, &jobs[ i ] );
        }
        for ( uint32_t i = 0; i < 4; ++i )
        {
            threads[ i ]->join();
            delete threads[ i ];
        }

        PrintEnergyBuffer( EAvgBuffer, s_RoughnessCount, 1, 1 );*/
    }
}


