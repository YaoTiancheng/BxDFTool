#pragma once

struct SRandomNumberGenerator;

struct SEnergyIntegral
{
    void Execute( uint32_t threadIndex, uint32_t localLaneIndex, uint32_t globalLaneIndex );

    SRandomNumberGenerator* m_Rngs;
    float    m_CosThetaInterval;
    float    m_AlphaInterval;
    float    m_EtaInterval;
    float    m_EtaBegin;
    float    m_EtaI;
    uint32_t m_SampleCount;
    uint32_t m_AlphaCount;
    float*   m_OutputBuffer;
};

struct SAverageEnergyIntegral
{
    void Execute( uint32_t threadIndex, uint32_t localLaneIndex, uint32_t globalLaneIndex );

    float*   m_EnergyBuffer;
    uint32_t m_CosThetaCount;
    float*   m_OutputBuffer;
};

struct SComputeInvCDF
{
    void Execute( uint32_t threadIndex, uint32_t localLaneIndex, uint32_t globalLaneIndex );

    float*   m_EnergyBuffer;
    uint32_t m_CosThetaCount;
    float*   m_OutputBuffer;
    float*   m_PDFScaleOutputBuffer;
    float    m_MaxPDFScale;
};