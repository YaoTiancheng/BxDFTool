#pragma once

#include <cstdint>

namespace NumericalIntegration
{
    float CompositeTrapezoidal( float a, float b, float* samples, uint32_t sampleCount );
}