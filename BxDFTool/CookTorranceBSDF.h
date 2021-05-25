#pragma once

#include <math.h>
#include "DirectXMath.h"
#include "LightingContext.h"

inline float EvaluateDielectricFresnel( float cosThetaI, float etaI, float etaT )
{
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

inline DirectX::XMFLOAT3 GGXSampleHemisphere( const DirectX::XMFLOAT2& sample, float alpha )
{
    float theta = atanf( alpha * sqrt( sample.x / ( 1.0f - sample.x ) ) );
    float phi = 2.0f * (float)M_PI * sample.y;

    float s = sin( theta );
    return DirectX::XMFLOAT3( cos( phi ) * s, sin( phi ) * s, cos( theta ) );
}

inline float EvaluateGGXMicrofacetDistribution( const DirectX::XMFLOAT3& m, float alpha )
{
    float alpha2 = alpha * alpha;
    float NdotM = m.z;
    float NdotM2 = NdotM * NdotM;
    float factor = NdotM2 * ( alpha2 - 1.0f ) + 1.0f;
    float denominator = factor * factor * (float)M_PI;
    return alpha2 / denominator;
}

inline float EvaluateGGXMicrofacetDistributionPdf( const DirectX::XMFLOAT3& m, float alpha )
{
    return EvaluateGGXMicrofacetDistribution( m, alpha ) * abs( m.z );
}

inline void SampleGGXMicrofacetDistribution( const DirectX::XMFLOAT2& sample, float alpha, DirectX::XMFLOAT3* m, float* pdf )
{
    *m = GGXSampleHemisphere( sample, alpha );
    *pdf = EvaluateGGXMicrofacetDistributionPdf( *m, alpha );
}

inline void SampleGGXMicrofacetDistribution( const DirectX::XMFLOAT2& sample, float alpha, DirectX::XMFLOAT3* m )
{
    *m = GGXSampleHemisphere( sample, alpha );
}

//
// GGX geometric shadowing
//

inline float EvaluateGGXGeometricShadowingOneDirection( float alpha2, const DirectX::XMFLOAT3& m, const DirectX::XMFLOAT3& w )
{
    // Ensure consistent orientation.
    // (Can't see backfacing microfacet normal from the front and vice versa)
    float WdotM;
    XMStoreFloat( &WdotM, XMVector3Dot( XMLoadFloat3( &w ), XMLoadFloat3( &m ) ) );
    if ( WdotM * w.z <= 0.0f )
        return 0.0f;

    float NdotW = abs( w.z );
    float denominator = sqrt( alpha2 + ( 1.0f - alpha2 ) * NdotW * NdotW ) + NdotW;
    return 2.0f * NdotW / denominator;
}

inline float EvaluateGGXGeometricShadowing( const DirectX::XMFLOAT3& wi, const DirectX::XMFLOAT3& wo, const DirectX::XMFLOAT3& m, float alpha )
{
    float alpha2 = alpha * alpha;
    return EvaluateGGXGeometricShadowingOneDirection( alpha2, m, wi ) * EvaluateGGXGeometricShadowingOneDirection( alpha2, m, wo );
}

inline float EvaluateSpecularBRDF( const DirectX::XMFLOAT3& wi, const DirectX::XMFLOAT3& wo )
{
    return 0.0f;
}

inline float EvaluateSpecularBRDFPdf( const DirectX::XMFLOAT3& wi, const DirectX::XMFLOAT3& wo )
{
    return 0.0f;
}

inline void SampleSpecularBRDF( const DirectX::XMFLOAT3& wo, float etaI, float etaT, DirectX::XMFLOAT3* wi, float* value, float* pdf, SLightingContext* lightingContext )
{
    *wi    = DirectX::XMFLOAT3( -wo.x, -wo.y, wo.z );
    *value = lightingContext->WOdotN > 0.0f ? EvaluateDielectricFresnel( lightingContext->WOdotN, etaI, etaT ) / wi->z : 0.0f;
    *pdf   = lightingContext->WOdotN > 0.0f ? 1.0f : 0.0f;

    lightingContext->WIdotN = wi->z;
    lightingContext->H = DirectX::XMFLOAT3( 0.0f, 0.0f, 1.0f );
    lightingContext->WOdotH = lightingContext->WOdotN;
}

inline float EvaluateSpecularBTDF( const DirectX::XMFLOAT3& wi, const DirectX::XMFLOAT3& wo )
{
    return 0.0f;
}

inline float EvaluateSpecularBTDFPdf( const DirectX::XMFLOAT3& wi, const DirectX::XMFLOAT3& wo )
{
    return 0.0f;
}

inline void SampleSpecularBTDF( const DirectX::XMFLOAT3& wo, float etaI, float etaT, DirectX::XMFLOAT3* wi, float* value, float* pdf, SLightingContext* lightingContext )
{
    *value = 0.0f;
    *pdf = 0.0f;

    DirectX::XMVECTOR xmWi = XMVector3Refract( XMVectorNegate( XMLoadFloat3( &wo ) ), DirectX::g_XMIdentityR2, etaI / etaT );
    XMStoreFloat3( wi, xmWi );
    if ( wi->x == 0.0f && wi->y == 0.0f && wi->z == 0.0f )
        return;

    lightingContext->WIdotN = wi->z;

    *value = lightingContext->WOdotN > 0.0f ? ( 1.0f - EvaluateDielectricFresnel( lightingContext->WIdotN, etaI, etaT ) ) : 0.0f;
    *value *= ( etaI * etaI ) / ( etaT * etaT );
    *value /= abs( lightingContext->WIdotN );

    *pdf = lightingContext->WOdotN > 0.0f ? 1.0f : 0.0f;
}

#define ALPHA_THRESHOLD 0.000196f

//
// Cook-Torrance microfacet BRDF
//

inline float EvaluateCookTorranceMircofacetBRDF( const DirectX::XMFLOAT3& wi, const DirectX::XMFLOAT3& wo, float alpha, float etaI, float etaT, const SLightingContext& lightingContext )
{
    if ( alpha >= ALPHA_THRESHOLD )
    {
        float WIdotN = lightingContext.WIdotN;
        float WOdotN = lightingContext.WOdotN;
        float WOdotM = lightingContext.WOdotH;
        if ( WIdotN <= 0.0f || WOdotN <= 0.0f || WOdotM <= 0.0f )
            return 0.0f;

        DirectX::XMFLOAT3 m = lightingContext.H;
        if ( m.x == 0.0f && m.y == 0.0f && m.z == 0.0f )
            return 0.0f;

        return EvaluateGGXMicrofacetDistribution( m, alpha ) * EvaluateGGXGeometricShadowing( wi, wo, m, alpha ) * EvaluateDielectricFresnel( std::min( 1.0f, WOdotM ), etaI, etaT ) / ( 4.0f * WIdotN * WOdotN );
    }
    else
    {
        return EvaluateSpecularBRDF( wi, wo );
    }
}

inline float EvaluateCookTorranceMicrofacetBRDFPdf( const DirectX::XMFLOAT3& wi, const DirectX::XMFLOAT3& wo, float alpha, const SLightingContext& lightingContext )
{
    if ( alpha >= ALPHA_THRESHOLD )
    {
        float WIdotN = lightingContext.WIdotN;
        if ( WIdotN <= 0.0f )
            return 0.0f;

        DirectX::XMFLOAT3 m = lightingContext.H;
        float WOdotM = lightingContext.WOdotH;
        float pdf = EvaluateGGXMicrofacetDistributionPdf( m, alpha );
        return pdf / ( 4.0f * WOdotM );
    }
    else
    {
        return EvaluateSpecularBRDFPdf( wi, wo );
    }
}

inline void SampleCookTorranceMicrofacetBRDF( const DirectX::XMFLOAT3& wo, const DirectX::XMFLOAT2& sample, float alpha, float etaI, float etaT, DirectX::XMFLOAT3* wi, float* value, float* pdf, bool* isDeltaBrdf, SLightingContext* lightingContext )
{
    if ( alpha >= ALPHA_THRESHOLD )
    {
        DirectX::XMFLOAT3 m;
        SampleGGXMicrofacetDistribution( sample, alpha, &m );
        XMStoreFloat3( wi, XMVectorNegate( XMVector3Reflect( XMLoadFloat3( &wo ), XMLoadFloat3( &m ) ) ) );

        lightingContext->H = m;
        XMStoreFloat( &lightingContext->WOdotH, XMVector3Dot( XMLoadFloat3( &wo ), XMLoadFloat3( &m ) ) );

        float WIdotN = wi->z;
        lightingContext->WIdotN = WIdotN;

        *value = EvaluateCookTorranceMircofacetBRDF( *wi, wo, alpha, etaI, etaT, *lightingContext );
        *pdf   = EvaluateCookTorranceMicrofacetBRDFPdf( *wi, wo, alpha, *lightingContext );
        *isDeltaBrdf = false;
    }
    else
    {
        SampleSpecularBRDF( wo, etaI, etaT, wi, value, pdf, lightingContext );
        *isDeltaBrdf = true;
    }
}

inline float EvaluateCookTorranceMircofacetBTDF( const DirectX::XMFLOAT3& wi, const DirectX::XMFLOAT3& wo, float alpha, float etaI, float etaT, const SLightingContext& lightingContext )
{
    if ( alpha >= ALPHA_THRESHOLD )
    {
        float WIdotN = lightingContext.WIdotN;
        float WOdotN = lightingContext.WOdotN;
        if ( WIdotN == 0.0f || WOdotN == 0.0f )
            return 0.0f;

        DirectX::XMFLOAT3 m;
        DirectX::XMVECTOR xmM = XMVectorAdd( XMVectorMultiply( XMLoadFloat3( &wo ), XMVectorReplicate( etaI ) )
                                           , XMVectorMultiply( XMLoadFloat3( &wi ), XMVectorReplicate( etaT ) ) );
        xmM = XMVector3Normalize( xmM );
        XMStoreFloat3( &m, xmM );
         
        if ( m.z < 0.0f ) m = XMFLOAT3( -m.x, -m.y, -m.z ); // Ensure same facing as the wo otherwise it will be rejected in G
        float WIdotM, WOdotM;
        XMStoreFloat( &WIdotM, XMVector3Dot( XMLoadFloat3( &wi ), xmM ) );
        XMStoreFloat( &WOdotM, XMVector3Dot( XMLoadFloat3( &wo ), xmM ) );
        float sqrtDenom = etaI * WOdotM + etaT * WIdotM;

        return ( 1.0f - EvaluateDielectricFresnel( WIdotM, etaI, etaT ) )
            * abs( EvaluateGGXMicrofacetDistribution( m, alpha ) * EvaluateGGXGeometricShadowing( wi, wo, m, alpha )
            * abs( WIdotM ) * abs( WOdotM ) 
            * etaI * etaI // etaT * etaT * ( ( etaI * etaI ) / ( etaT * etaT ) )
            / ( WOdotN * WIdotN * sqrtDenom * sqrtDenom ) );
    }
    else
    {
        return EvaluateSpecularBTDF( wi, wo );
    }
}

inline float EvaluateCookTorranceMicrofacetBTDFPdf( const DirectX::XMFLOAT3& wi, const DirectX::XMFLOAT3& wo, float alpha, float etaI, float etaT, const SLightingContext& lightingContext )
{
    if ( alpha >= ALPHA_THRESHOLD )
    {
        DirectX::XMFLOAT3 m;
        DirectX::XMVECTOR xmM = XMVectorAdd( XMVectorMultiply( XMLoadFloat3( &wo ), XMVectorReplicate( etaI ) )
                                           , XMVectorMultiply( XMLoadFloat3( &wi ), XMVectorReplicate( etaT ) ) );
        xmM = XMVector3Normalize( xmM );
        XMStoreFloat3( &m, xmM );

        if ( m.z < 0.0f ) m = XMFLOAT3( -m.x, -m.y, -m.z ); // Ensure same facing as the wo otherwise it will be rejected in G
        float WIdotM, WOdotM;
        XMStoreFloat( &WIdotM, XMVector3Dot( XMLoadFloat3( &wi ), xmM ) );
        XMStoreFloat( &WOdotM, XMVector3Dot( XMLoadFloat3( &wo ), xmM ) );
        float sqrtDenom = etaI * WOdotM + etaT * WIdotM;

        float dwh_dwi = abs( ( etaT * etaT * WIdotM ) / ( sqrtDenom * sqrtDenom ) );
        return EvaluateGGXMicrofacetDistributionPdf( m, alpha ) * dwh_dwi;
    }
    else
    {
        return EvaluateSpecularBTDFPdf( wi, wo );
    }
}

inline void SampleCookTorranceMicrofacetBTDF( const DirectX::XMFLOAT3& wo, const DirectX::XMFLOAT2& sample, float alpha, float etaI, float etaT, DirectX::XMFLOAT3* wi, float* value, float* pdf, bool* isDeltaBtdf, SLightingContext* lightingContext )
{
    *value = 0.0f;
    *wi    = XMFLOAT3( 0.0f, 0.0f, 0.0f );
    *pdf   = 0.0f;

    if ( alpha >= ALPHA_THRESHOLD )
    {
        float WOdotN = lightingContext->WOdotN;
        if ( WOdotN == 0.0f )
            return;

        DirectX::XMFLOAT3 m;
        SampleGGXMicrofacetDistribution( sample, alpha, &m );
        DirectX::XMVECTOR xmWi = XMVector3Refract( XMVectorNegate( XMLoadFloat3( &wo ) ), XMLoadFloat3( &m ), etaI / etaT );
        XMStoreFloat3( wi, xmWi );
        if ( wi->x == 0.0f && wi->y == 0.0f && wi->z == 0.0f )
            return;

        lightingContext->H = m;
        XMStoreFloat( &lightingContext->WOdotH, XMVector3Dot( XMLoadFloat3( &wo ), XMLoadFloat3( &m ) ) );

        float WIdotN = wi->z;
        lightingContext->WIdotN = WIdotN;

        *value = EvaluateCookTorranceMircofacetBTDF( *wi, wo, alpha, etaI, etaT, *lightingContext );
        *pdf   = EvaluateCookTorranceMicrofacetBTDFPdf( *wi, wo, alpha, etaI, etaT, *lightingContext );

        *isDeltaBtdf = false;
    }
    else
    {
        SampleSpecularBTDF( wo, etaI, etaT, wi, value, pdf, lightingContext );
        *isDeltaBtdf = true;
    }
}

inline float EvaluateCookTorranceMicrofacetBSDF( const DirectX::XMFLOAT3& wi, const DirectX::XMFLOAT3& wo, float alpha, float etaI, float etaT, const SLightingContext& lightingContext )
{
    if ( wi.z < 0.0f )
        return EvaluateCookTorranceMircofacetBTDF( wi, wo, alpha, etaI, etaT, lightingContext );
    else
        return EvaluateCookTorranceMircofacetBRDF( wi, wo, alpha, etaI, etaT, lightingContext );
}

inline float EvaluateCookTorranceMicrofacetBSDFPdf( const DirectX::XMFLOAT3& wi, const DirectX::XMFLOAT3& wo, float alpha, float etaI, float etaT, const SLightingContext& lightingContext )
{
    if ( wi.z < 0.0f )
        return EvaluateCookTorranceMicrofacetBTDFPdf( wi, wo, alpha, etaI, etaT, lightingContext );
    else
        return EvaluateCookTorranceMicrofacetBRDFPdf( wi, wo, alpha, lightingContext );
}

inline void SampleCookTorranceMicrofacetBSDF( const DirectX::XMFLOAT3& wo, float selectionSample, const DirectX::XMFLOAT2& bxdfSample, float alpha, float etaI, float etaT, DirectX::XMFLOAT3* wi, float* value, float* pdf, bool* isDeltaBxdf, SLightingContext* lightingContext )
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
    static void Sample( const DirectX::XMFLOAT3& wo, float selectionSample, const DirectX::XMFLOAT2& bxdfSample, float alpha, float etaI, float etaT, DirectX::XMFLOAT3* wi, float* value, float* pdf, bool* isDeltaBxdf, SLightingContext* lightingContext )
    {
        SampleCookTorranceMicrofacetBRDF( wo, bxdfSample, alpha, etaI, etaT, wi, value, pdf, isDeltaBxdf, lightingContext );
    }

    static const bool s_NeedSelectionSample = false;
};


class CookTorranceMicrofacetBSDF
{
public:
    static void Sample( const DirectX::XMFLOAT3& wo, float selectionSample, const DirectX::XMFLOAT2& bxdfSample, float alpha, float etaI, float etaT, DirectX::XMFLOAT3* wi, float* value, float* pdf, bool* isDeltaBxdf, SLightingContext* lightingContext )
    {
        SampleCookTorranceMicrofacetBSDF( wo, selectionSample, bxdfSample, alpha, etaI, etaT, wi, value, pdf, isDeltaBxdf, lightingContext );
    }

    static const bool s_NeedSelectionSample = true;
};