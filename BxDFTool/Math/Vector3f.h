#pragma once

struct Vector3f
{
    Vector3f() = default;
    Vector3f( float x, float y, float z ) : x( x ), y( y ), z( z ) {}
    Vector3f( const Vector3f& other ) : x( other.x ), y( other.y ), z( other.z ) {}
    Vector3f( const Vector3f&& other ) noexcept : x( other.x ), y( other.y ), z( other.z ) {}

    Vector3f operator+( const Vector3f& other ) const
    {
        return Vector3f( x + other.x, y + other.y, z + other.z );
    }

    Vector3f operator-( const Vector3f& other ) const
    {
        return Vector3f( x - other.x, y - other.y, z - other.z );
    }

    Vector3f operator*( const Vector3f& other ) const
    {
        return Vector3f( x * other.x, y * other.y, z * other.z );
    }

    Vector3f operator/( const Vector3f& other ) const
    {
        return Vector3f( x / other.x, y / other.y, z / other.z );
    }

    Vector3f operator+( float s ) const
    {
        return Vector3f( x + s, y + s, z + s );
    }

    Vector3f operator-( float s ) const
    {
        return Vector3f( x - s, y - s, z - s );
    }

    Vector3f operator*( float s ) const
    {
        return Vector3f( x * s, y * s, z * s );
    }

    Vector3f operator/( float s ) const
    {
        return Vector3f( x / s, y / s, z / s );
    }

    Vector3f& operator=( const Vector3f& other )
    {
        x = other.x;
        y = other.y;
        z = other.z;
        return *this;
    }

    static float Dot( const Vector3f& lhs, const Vector3f& rhs )
    {
        return lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z;
    }

    static Vector3f Cross( const Vector3f& lhs, const Vector3f& rhs )
    {
        return Vector3f( lhs.y * rhs.z - lhs.z * rhs.y
                       , lhs.z * rhs.x - lhs.x * rhs.z
                       , lhs.x * rhs.y - lhs.y * rhs.x );
    }

    union
    {
        struct { float x; float y; float z; };
        struct { float m[ 3 ]; };
    };
};

inline Vector3f operator-( const Vector3f& v )
{
    return Vector3f( -v.x, -v.y, -v.z );
}