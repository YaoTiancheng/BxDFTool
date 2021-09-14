#pragma once

#include <DirectXMath.h>
#include "Vector2f.h"
#include "Vector3f.h"

namespace MathHelpers
{
    inline DirectX::XMVECTOR XM_CALLCONV Vec3fToVector( const Vector3f& v )
    {
        return XMLoadFloat3( (DirectX::XMFLOAT3*)&v );
    }

    inline Vector3f XM_CALLCONV VectorToVec3f( FXMVECTOR v )
    {
        Vector3f result;
        XMStoreFloat3( (DirectX::XMFLOAT3*)&result, v );
        return result;
    }

    inline Vector3f Reflect( const Vector3f& i, const Vector3f& n )
    {
        return VectorToVec3f( DirectX::XMVector3Reflect( Vec3fToVector( i ), Vec3fToVector( n ) ) );
    }

    inline Vector3f Refract( const Vector3f& i, const Vector3f& n, float refractionIndex )
    {
        return VectorToVec3f( DirectX::XMVector3Refract( Vec3fToVector( i ), Vec3fToVector( n ), refractionIndex ) );
    }

    inline float Clamp( float v, float min, float max )
    {
        return std::max( std::min( v, max ), min );
    }
}

using float2 = Vector2f;
using float3 = Vector3f;