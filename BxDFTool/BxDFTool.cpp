// BRDFIntegral.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "PCH.h"
#include "RandomNumberGenerator.h"
#include "EnergyIntegral.h"
#include "MultiProcessing/MPDispatcher.h"

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
                value = std::roundf( std::min( 1.0f, value ) * 65535.0f );
                //printf( "%9.7ff, ", *pBuffer );
                printf( "%5d, ", (uint16_t)value );
                pBuffer++;
            }
            printf( "\n" );
        }
        printf( "\n" );
    }
}

int main()
{
    uint32_t alphaCount = 16;
    uint32_t cosThetaCount = 32;
    uint32_t etaCount = 16;
    uint32_t sampleCount = 960000;
    bool     outputEnergyBuffer = true;
    bool     outputAverageEnergyBuffer = false;
    bool     outputInverseCDF = true;
    float    etaI = 1.0f;
    float    etaTBegin = 1.0f;
    float    etaTEnd = 3.0f;
    bool     invert = true;

    float* energyBuffer = nullptr;
    {
        uint32_t threadSize = cosThetaCount;
        uint32_t threadCount = alphaCount * etaCount;
        energyBuffer = new float[ threadSize * threadCount ];
        SRandomNumberGenerator* rngs = new SRandomNumberGenerator[ threadCount ];

        SEnergyIntegral integral;
        integral.m_AlphaInterval    = 1.0f / ( alphaCount - 1 );
        integral.m_CosThetaInterval = 1.0f / ( cosThetaCount - 1 );
        integral.m_EtaBegin         = etaTBegin;
        integral.m_EtaInterval      = etaCount > 1 ? ( etaTEnd - integral.m_EtaBegin ) / ( etaCount - 1 ) : 0.0f;
        integral.m_EtaI             = etaI;
        integral.m_AlphaCount       = alphaCount;
        integral.m_Invert           = invert;
        integral.m_OutputBuffer     = energyBuffer;
        integral.m_Rngs             = rngs;
        integral.m_SampleCount      = sampleCount;

        using std::placeholders::_1;
        using std::placeholders::_2;
        using std::placeholders::_3;
        MP::LaneFunctionType laneFunction = std::bind( &SEnergyIntegral::Execute_CookTorranceMicrofacetBSDF, &integral, _1, _2, _3 );
        MP::Dispatch( threadSize, threadCount, laneFunction );

        delete[] rngs;

        if ( outputEnergyBuffer )
        {
            PrintEnergyBuffer( energyBuffer, cosThetaCount, alphaCount, etaCount );
        }
    }

    if ( outputAverageEnergyBuffer )
    {
        uint32_t threadSize = 1;
        uint32_t threadCount = alphaCount * etaCount;

        float* outputBuffer = new float[ threadSize * threadCount ];
        SAverageEnergyIntegral integral;
        integral.m_EnergyBuffer  = energyBuffer;
        integral.m_CosThetaCount = cosThetaCount;
        integral.m_OutputBuffer  = outputBuffer;

        using std::placeholders::_1;
        using std::placeholders::_2;
        using std::placeholders::_3;
        MP::LaneFunctionType laneFunction = std::bind( &SAverageEnergyIntegral::Execute, &integral, _1, _2, _3 );
        MP::Dispatch( threadSize, threadCount, laneFunction );

        PrintEnergyBuffer( outputBuffer, alphaCount, etaCount, 1 );

        delete[] outputBuffer;
    }

    if ( outputInverseCDF )
    {
        uint32_t threadSize = 1;
        uint32_t threadCount = alphaCount * etaCount;

        float* outputBuffer = new float[ cosThetaCount * alphaCount * etaCount ];
        float* PDFScaleOutputBuffer = new float[ alphaCount * etaCount ];
        SComputeInvCDF invCDF;
        invCDF.m_CosThetaCount = cosThetaCount;
        invCDF.m_EnergyBuffer  = energyBuffer;
        invCDF.m_MaxPDFScale   = 2.0f;
        invCDF.m_OutputBuffer  = outputBuffer;
        invCDF.m_PDFScaleOutputBuffer = PDFScaleOutputBuffer;

        using std::placeholders::_1;
        using std::placeholders::_2;
        using std::placeholders::_3;
        MP::LaneFunctionType laneFunction = std::bind( &SComputeInvCDF::Execute, &invCDF, _1, _2, _3 );
        MP::Dispatch( threadSize, threadCount, laneFunction );

        PrintEnergyBuffer( outputBuffer, cosThetaCount, alphaCount, etaCount );
        PrintEnergyBuffer( PDFScaleOutputBuffer, alphaCount, etaCount, 1 );

        delete[] outputBuffer;
        delete[] PDFScaleOutputBuffer;
    }

    delete[] energyBuffer;
}


