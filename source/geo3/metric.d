module geo3.metric;

import geo3.point :
    Point3;

import geo3.scalar :
    MetricScalar,
    isGeoScalar;

import std.math.algebraic :
    hypot;

MetricScalar!T distance(T)(
    Point3!T a,
    Point3!T b
)
    pure nothrow @safe @nogc
if (isGeoScalar!T)
{
    alias M = MetricScalar!T;

    const M dx =
        cast(M) a.x -
        cast(M) b.x;

    const M dy =
        cast(M) a.y -
        cast(M) b.y;

    const M dz =
        cast(M) a.z -
        cast(M) b.z;

    return hypot(
        hypot(
            dx,
            dy
        ),
        dz
    );
}

@safe unittest
{
    const a =
        Point3!double(
            0.0,
            0.0,
            0.0
        );

    const b =
        Point3!double(
            2.0,
            3.0,
            6.0
        );

    assert(
        distance(a, b) ==
        7.0
    );
}
