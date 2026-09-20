module geo3.point;

import geo3.scalar :
    isGeoScalar;

struct Point3(T)
if (isGeoScalar!T)
{
    private T _x = T(0);
    private T _y = T(0);
    private T _z = T(0);

    this(
        T x,
        T y,
        T z
    )
        pure nothrow @safe @nogc
    {
        _x = x;
        _y = y;
        _z = z;
    }

    @property T x() const
        pure nothrow @safe @nogc
    {
        return _x;
    }

    @property T y() const
        pure nothrow @safe @nogc
    {
        return _y;
    }

    @property T z() const
        pure nothrow @safe @nogc
    {
        return _z;
    }
}

@safe unittest
{
    const p =
        Point3!double(
            1.0,
            2.0,
            3.0
        );

    assert(p.x == 1.0);
    assert(p.y == 2.0);
    assert(p.z == 3.0);
}
