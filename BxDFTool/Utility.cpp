#include "Utility.h"

float NumericalIntegration::CompositeTrapezoidal( float a, float b, float* samples, uint32_t sampleCount )
{
    uint32_t n = sampleCount - 1;
    float fa = samples[ 0 ];
    float fb = samples[ sampleCount - 1 ];
    float sum = 0;
    for ( uint32_t i = 1; i < sampleCount - 1; ++i )
    {
        sum += samples[ i ];
    }
    return ( sum + ( fa + fb ) * 0.5f ) * ( b - a ) / n;;
}