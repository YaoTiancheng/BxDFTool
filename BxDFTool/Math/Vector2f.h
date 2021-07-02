#pragma once

struct Vector2f
{
    Vector2f() = default;
    Vector2f( float x, float y ) : x( x ), y( y ) {}
    Vector2f( const Vector2f& other ) : x( other.x ), y( other.y ) {}
    Vector2f( const Vector2f&& other ) noexcept : x( other.x ), y( other.y ) {}

    Vector2f operator+( const Vector2f& other ) const
    {
        return Vector2f( x + other.x, y + other.y );
    }

    Vector2f operator-( const Vector2f& other ) const
    {
        return Vector2f( x - other.x, y - other.y );
    }

    Vector2f operator*( const Vector2f& other ) const
    {
        return Vector2f( x * other.x, y * other.y );
    }

    Vector2f operator/( const Vector2f& other ) const
    {
        return Vector2f( x / other.x, y / other.y );
    }

    Vector2f operator+( float s ) const
    {
        return Vector2f( x + s, y + s );
    }

    Vector2f operator-( float s ) const
    {
        return Vector2f( x - s, y - s );
    }

    Vector2f operator*( float s ) const
    {
        return Vector2f( x * s, y * s );
    }

    Vector2f operator/( float s ) const
    {
        return Vector2f( x / s, y / s );
    }

    Vector2f& operator=( const Vector2f& other )
    {
        x = other.x;
        y = other.y;
        return *this;
    }

    union
    {
        struct { float x; float y; };
        struct { float m[ 2 ]; };
    };
};