#pragma once

#include <math.h>
#include "Math/MathHelpers.h"
#include "LightingContext.h"

inline float EvaluateDielectricFresnel( float cosThetaI, float etaI, float etaT )
{
    cosThetaI = std::min( 1.0f, std::max( -1.0f, cosThetaI ) );

    if ( cosThetaI < 0.0f )
    {
        float etaTemp = etaI;
        etaI = etaT;
        etaT = etaTemp;
        cosThetaI = -cosThetaI;
    }

    float sinThetaI = sqrtf( 1.0f - cosThetaI * cosThetaI );
    float sinThetaT = etaI / etaT * sinThetaI;
    if ( sinThetaT >= 1.0f )
    {
        return 1.0f;
    }
    float cosThetaT = sqrtf( 1.0f - sinThetaT * sinThetaT );
    float Rparl = ( ( etaT * cosThetaI ) - ( etaI * cosThetaT ) ) /
        ( ( etaT * cosThetaI ) + ( etaI * cosThetaT ) );
    float Rperp = ( ( etaI * cosThetaI ) - ( etaT * cosThetaT ) ) /
        ( ( etaI * cosThetaI ) + ( etaT * cosThetaT ) );
    return ( Rparl * Rparl + Rperp * Rperp ) * 0.5f;
}

inline float3 GGXSampleHemisphere( const float2& sample, float alpha )
{
    float theta = atanf( alpha * sqrt( sample.x / ( 1.0f - sample.x ) ) );
    float phi = 2.0f * (float)M_PI * sample.y;

    float s = sin( theta );
    return float3( cos( phi ) * s, sin( phi ) * s, cos( theta ) );
}

inline float EvaluateGGXMicrofacetDistribution( const float3& m, float alpha )
{
    float alpha2 = alpha * alpha;
    float NdotM = m.z;
    float NdotM2 = NdotM * NdotM;
    float factor = NdotM2 * ( alpha2 - 1.0f ) + 1.0f;
    float denominator = factor * factor * (float)M_PI;
    return alpha2 / denominator;
}

inline float EvaluateGGXMicrofacetDistributionPdf( const float3& m, float alpha )
{
    return EvaluateGGXMicrofacetDistribution( m, alpha ) * abs( m.z );
}

inline void SampleGGXMicrofacetDistribution( const float2& sample, float alpha, float3* m, float* pdf )
{
    *m = GGXSampleHemisphere( sample, alpha );
    *pdf = EvaluateGGXMicrofacetDistributionPdf( *m, alpha );
}

inline void SampleGGXMicrofacetDistribution( const float2& sample, float alpha, float3* m )
{
    *m = GGXSampleHemisphere( sample, alpha );
}

//
// GGX geometric shadowing
//

inline float EvaluateGGXGeometricShadowingOneDirection( float alpha2, const float3& m, const float3& w )
{
    // Ensure consistent orientation.
    // (Can't see backfacing microfacet normal from the front and vice versa)
    float WdotM = Vector3f::Dot( w, m );
    if ( WdotM * w.z <= 0.0f )
        return 0.0f;

    float NdotW = abs( w.z );
    float denominator = sqrt( alpha2 + ( 1.0f - alpha2 ) * NdotW * NdotW ) + NdotW;
    return 2.0f * NdotW / denominator;
}

inline float EvaluateGGXGeometricShadowing( const float3& wi, const float3& wo, const float3& m, float alpha )
{
    float alpha2 = alpha * alpha;
    return EvaluateGGXGeometricShadowingOneDirection( alpha2, m, wi ) * EvaluateGGXGeometricShadowingOneDirection( alpha2, m, wo );
}

#define ALPHA_THRESHOLD 0.000196f

//
// Cook-Torrance microfacet BRDF
//

inline void SampleCookTorranceMicrofacetBRDF( const float3& wo, const float2& sample, float alpha, float etaI, float etaT, float3* wi, float* value, float* pdf, bool* isDeltaBrdf, SLightingContext* lightingContext )
{
    *value = 0.0f;
    *wi    = float3( 0.0f, 0.0f, 0.0f );
    *pdf   = 0.0f;
    *isDeltaBrdf = false;

    if ( wo.z <= 0.0f )
        return;

    if ( alpha >= ALPHA_THRESHOLD )
    {
        float3 m;
        SampleGGXMicrofacetDistribution( sample, alpha, &m );
        *wi = -MathHelpers::Reflect( wo, m );

        lightingContext->H = m;
        lightingContext->WOdotH = Vector3f::Dot( wo, m );

        float WIdotN = wi->z;
        lightingContext->WIdotN = WIdotN;

        float WOdotN = lightingContext->WOdotN;
        float WOdotM = lightingContext->WOdotH;
        if ( WIdotN <= 0.0f || WOdotM <= 0.0f )
            return;

        *value = EvaluateGGXMicrofacetDistribution( m, alpha ) * EvaluateGGXGeometricShadowing( *wi, wo, m, alpha ) * EvaluateDielectricFresnel( std::min( 1.0f, WOdotM ), etaI, etaT ) / ( 4.0f * WIdotN * WOdotN );
        *pdf = EvaluateGGXMicrofacetDistributionPdf( m, alpha ) / ( 4.0f * WOdotM );
    }
    else
    {
        *wi = float3( -wo.x, -wo.y, wo.z );

        lightingContext->WIdotN = wi->z;
        lightingContext->H = float3( 0.0f, 0.0f, 1.0f );
        lightingContext->WOdotH = lightingContext->WOdotN;

        if ( wi->z == 0.0f )
            return;

        *value = EvaluateDielectricFresnel( lightingContext->WOdotN, etaI, etaT ) / wi->z;
        *pdf = 1.0f;
        *isDeltaBrdf = true;
    }
}

inline void SampleCookTorranceMicrofacetBTDF( const float3& wo, const float2& sample, float alpha, float etaI, float etaT, float3* wi, float* value, float* pdf, bool* isDeltaBtdf, SLightingContext* lightingContext )
{
    *value = 0.0f;
    *wi    = float3( 0.0f, 0.0f, 0.0f );
    *pdf   = 0.0f;
    *isDeltaBtdf = false;

    if ( wo.z <= 0.0f )
        return;

    if ( alpha >= ALPHA_THRESHOLD )
    {
        float WOdotN = lightingContext->WOdotN;

        float3 m;
        SampleGGXMicrofacetDistribution( sample, alpha, &m );

        *wi = MathHelpers::Refract( -wo, m, etaI / etaT );

        lightingContext->H = m;
        lightingContext->WOdotH = Vector3f::Dot( wo, m );

        float WIdotN = wi->z;
        lightingContext->WIdotN = WIdotN;

        if ( WIdotN == 0.0f )
            return;

        float WIdotM = Vector3f::Dot( *wi, m );
        float WOdotM = lightingContext->WOdotH;
        float sqrtDenom = etaI * WOdotM + etaT * WIdotM;

        *value = ( 1.0f - EvaluateDielectricFresnel( WIdotM, etaI, etaT ) )
            * abs( EvaluateGGXMicrofacetDistribution( m, alpha ) * EvaluateGGXGeometricShadowing( *wi, wo, m, alpha )
            * abs( WIdotM ) * abs( WOdotM )
            //* etaI * etaI // etaT * etaT * ( ( etaI * etaI ) / ( etaT * etaT ) )
            * etaT * etaT
            / ( WOdotN * WIdotN * sqrtDenom * sqrtDenom ) );
        assert( *value >= 0.0f );

        float dwh_dwi = abs( ( etaT * etaT * WIdotM ) / ( sqrtDenom * sqrtDenom ) );
        *pdf = EvaluateGGXMicrofacetDistributionPdf( m, alpha ) * dwh_dwi;
    }
    else
    {
        *wi = MathHelpers::Refract( -wo, float3( 0.0f, 0.0f, 1.0f ), etaI / etaT );

        lightingContext->WIdotN = wi->z;
        if ( wi->z == 0.0f )
            return;

        *value = 1.0f - EvaluateDielectricFresnel( lightingContext->WIdotN, etaI, etaT );
        //*value *= ( etaI * etaI ) / ( etaT * etaT );
        *value /= abs( lightingContext->WIdotN );
        *pdf = 1.0f;
        *isDeltaBtdf = true;
    }
}

inline void SampleCookTorranceMicrofacetBSDF( const float3& wo, float selectionSample, const float2& bxdfSample, float alpha, float etaI, float etaT, float3* wi, float* value, float* pdf, bool* isDeltaBxdf, SLightingContext* lightingContext )
{
    if ( selectionSample < 0.5f )
        SampleCookTorranceMicrofacetBRDF( wo, bxdfSample, alpha, etaI, etaT, wi, value, pdf, isDeltaBxdf, lightingContext );
    else
        SampleCookTorranceMicrofacetBTDF( wo, bxdfSample, alpha, etaI, etaT, wi, value, pdf, isDeltaBxdf, lightingContext );

    *pdf *= 0.5f;
}


class CookTorranceMicrofacetBRDF
{
public:
    static void Sample( const float3& wo, float selectionSample, const float2& bxdfSample, float alpha, float etaI, float etaT, float3* wi, float* value, float* pdf, bool* isDeltaBxdf, SLightingContext* lightingContext )
    {
        SampleCookTorranceMicrofacetBRDF( wo, bxdfSample, alpha, etaI, etaT, wi, value, pdf, isDeltaBxdf, lightingContext );
    }

    static const bool s_NeedSelectionSample = false;
};


class CookTorranceMicrofacetBTDF
{
public:
    static void Sample( const float3& wo, float selectionSample, const float2& bxdfSample, float alpha, float etaI, float etaT, float3* wi, float* value, float* pdf, bool* isDeltaBxdf, SLightingContext* lightingContext )
    {
        SampleCookTorranceMicrofacetBTDF( wo, bxdfSample, alpha, etaI, etaT, wi, value, pdf, isDeltaBxdf, lightingContext );
    }

    static const bool s_NeedSelectionSample = false;
};


class CookTorranceMicrofacetBSDF
{
public:
    static void Sample( const float3& wo, float selectionSample, const float2& bxdfSample, float alpha, float etaI, float etaT, float3* wi, float* value, float* pdf, bool* isDeltaBxdf, SLightingContext* lightingContext )
    {
        SampleCookTorranceMicrofacetBSDF( wo, selectionSample, bxdfSample, alpha, etaI, etaT, wi, value, pdf, isDeltaBxdf, lightingContext );
    }

    static const bool s_NeedSelectionSample = true;
};