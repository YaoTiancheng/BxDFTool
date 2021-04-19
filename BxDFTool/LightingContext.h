#pragma once

#include "DirectXMath.h"

struct SLightingContext
{
    void Init( const DirectX::XMFLOAT3& wo )
    {
        memset( this, 0, sizeof( SLightingContext ) );
        WOdotN = wo.z;
    }

    float WOdotN;
    float WIdotN;
    DirectX::XMFLOAT3 H;
    float WOdotH;
};