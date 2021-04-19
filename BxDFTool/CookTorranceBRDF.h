#pragma once

#include <math.h>
#include <algorithm>
#define _USE_MATH_DEFINES
#include "DirectXMath.h"
#include "LightingContext.h"

DirectX::XMFLOAT3 GGXSampleHemisphere( const DirectX::XMFLOAT2& sample, float alpha )
{
    float theta = atan( alpha * sqrt( sample.x / ( 1.0f - sample.x ) ) );
    float phi = 2.0f * (float)M_PI * sample.y;

    float s = sin( theta );
    return DirectX::XMFLOAT3( cos( phi ) * s, sin( phi ) * s, cos( theta ) );
}

float EvaluateGGXMicrofacetDistribution( const DirectX::XMFLOAT3& m, float alpha )
{
    float alpha2 = alpha * alpha;
    float NdotM = m.z;
    float NdotM2 = NdotM * NdotM;
    float factor = NdotM2 * ( alpha2 - 1.0f ) + 1.0f;
    float denominator = factor * factor * (float)M_PI;
    return alpha2 / denominator;
}

float EvaluateGGXMicrofacetDistributionPdf( const DirectX::XMFLOAT3& m, float alpha )
{
    return EvaluateGGXMicrofacetDistribution( m, alpha ) * m.z;
}

void SampleGGXMicrofacetDistribution( const DirectX::XMFLOAT2& sample, float alpha, DirectX::XMFLOAT3* m, float* pdf )
{
    *m = GGXSampleHemisphere( sample, alpha );
    *pdf = EvaluateGGXMicrofacetDistributionPdf( *m, alpha );
}

void SampleGGXMicrofacetDistribution( const DirectX::XMFLOAT2& sample, float alpha, DirectX::XMFLOAT3* m )
{
    *m = GGXSampleHemisphere( sample, alpha );
}

//
// GGX geometric shadowing
//

float EvaluateGGXGeometricShadowingOneDirection( float alpha2, const DirectX::XMFLOAT3& w )
{
    float NdotW = abs( w.z );
    float denominator = sqrtf( alpha2 + ( 1.0f - alpha2 ) * NdotW * NdotW ) + NdotW;
    return 2.0f * NdotW / denominator;
}

float EvaluateGGXGeometricShadowing( const DirectX::XMFLOAT3& wi, const DirectX::XMFLOAT3& wo, float alpha )
{
    float alpha2 = alpha * alpha;
    return EvaluateGGXGeometricShadowingOneDirection( alpha2, wi ) * EvaluateGGXGeometricShadowingOneDirection( alpha2, wo );
}

float EvaluateDielectricFresnel( float WOdotM, float etaI, float etaT )
{
    float cosThetaI = WOdotM;
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

//
// Cook-Torrance microfacet BRDF
//

float EvaluateCookTorranceMircofacetBRDF( const DirectX::XMFLOAT3& wi, const DirectX::XMFLOAT3& wo, float alpha, float etaI, float etaT, const SLightingContext& lightingContext )
{
    float WIdotN = lightingContext.WIdotN;
    float WOdotN = lightingContext.WOdotN;
    float WOdotM = lightingContext.WOdotH;
    if ( WIdotN <= 0.0f || WOdotN <= 0.0f || WOdotM <= 0.0f )
        return 0.0f;

    const DirectX::XMFLOAT3& m = lightingContext.H;
    if ( m.x == 0.0f && m.y == 0.0f && m.z == 0.0f )
        return 0.0f;

    return EvaluateGGXMicrofacetDistribution( m, alpha ) * EvaluateGGXGeometricShadowing( wi, wo, alpha ) * EvaluateDielectricFresnel( std::min( 1.0f, WOdotM ), etaI, etaT ) / ( 4.0f * WIdotN * WOdotN );
}

float EvaluateCookTorranceMicrofacetBRDFPdf( const DirectX::XMFLOAT3& wi, const DirectX::XMFLOAT3& wo, float alpha, const SLightingContext& lightingContext )
{
    float WIdotN = lightingContext.WIdotN;
    if ( WIdotN <= 0.0f )
        return 0.0f;

    const DirectX::XMFLOAT3& m = lightingContext.H;
    float pdf = EvaluateGGXMicrofacetDistributionPdf( m, alpha );
    float WOdotM = lightingContext.WOdotH;
    return pdf / ( 4.0f * WOdotM );
}

void SampleCookTorranceMicrofacetBRDF( const DirectX::XMFLOAT3& wo, const DirectX::XMFLOAT2& sample, float alpha, float etaI, float etaT, DirectX::XMFLOAT3* wi, float* value, float* pdf, SLightingContext* lightingContext )
{
    if ( alpha >= 0.001f )
    {
        DirectX::XMFLOAT3 m;
        SampleGGXMicrofacetDistribution( sample, alpha, &m );
        XMStoreFloat3( wi, DirectX::XMVectorNegate( DirectX::XMVector3Reflect( DirectX::XMLoadFloat3( &wo ), DirectX::XMLoadFloat3( &m ) ) ) );

        lightingContext->H = m;
        DirectX::XMStoreFloat( &lightingContext->WOdotH, DirectX::XMVector3Dot( DirectX::XMLoadFloat3( &wo ), DirectX::XMLoadFloat3( &m ) ) );

        float WIdotN = wi->z;
        lightingContext->WIdotN = WIdotN;

        *value = EvaluateCookTorranceMircofacetBRDF( *wi, wo, alpha, etaI, etaT, *lightingContext );
        *pdf = EvaluateCookTorranceMicrofacetBRDFPdf( *wi, wo, alpha, *lightingContext );
    }
    else
    {
        *wi = DirectX::XMFLOAT3( -wo.x, -wo.y, wo.z );
        *value = EvaluateDielectricFresnel( lightingContext->WOdotN, etaI, etaT ) / wi->z;
        *pdf = 1.0f;
    }
}

