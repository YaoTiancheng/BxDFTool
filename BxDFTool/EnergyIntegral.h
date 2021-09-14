#pragma once

struct SRandomNumberGenerator;

struct SEnergyIntegral
{
private:
    template <typename BxDF>
    void Execute( uint32_t threadIndex, uint32_t localLaneIndex, uint32_t globalLaneIndex );


public:
    void Execute_CookTorranceMicrofacetBRDF( uint32_t threadIndex, uint32_t localLaneIndex, uint32_t globalLaneIndex );
    void Execute_CookTorranceMicrofacetBTDF( uint32_t threadIndex, uint32_t localLaneIndex, uint32_t globalLaneIndex );
    void Execute_CookTorranceMicrofacetBSDF( uint32_t threadIndex, uint32_t localLaneIndex, uint32_t globalLaneIndex );
    void Execute_FresnelConductor( uint32_t threadIndex, uint32_t localLaneIndex, uint32_t globalLaneIndex );

    SRandomNumberGenerator* m_Rngs;
    float    m_CosThetaInterval;
    float    m_AlphaInterval;
    float    m_EtaInterval;
    float    m_EtaBegin;
    float    m_EtaI;
    float    m_kInterval;
    uint32_t m_SampleCount;
    uint32_t m_AlphaCount;
    bool     m_Invert;
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