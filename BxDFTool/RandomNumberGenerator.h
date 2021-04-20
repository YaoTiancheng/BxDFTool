#pragma once

#include <random>

struct SRandomNumberGenerator
{
    SRandomNumberGenerator()
        : m_Distribution( 0.0f, 1.0f )
    {
    }

    float Next()
    {
        return m_Distribution( m_Generator );
    }

    std::default_random_engine            m_Generator;
    std::uniform_real_distribution<float> m_Distribution;
};