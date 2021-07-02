#pragma once

#include "Math/Vector3f.h"

struct SLightingContext
{
    void Init( const float3& wo )
    {
        memset( this, 0, sizeof( SLightingContext ) );
        WOdotN = wo.z;
    }

    float WOdotN;
    float WIdotN;
    float3 H;
    float WOdotH;
};