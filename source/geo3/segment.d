/**
 * Three-dimensional closed line-segment primitives.
 *
 * Authors:
 *     Alexander Bernardi
 *
 * Copyright:
 *     Copyright © 2026 Alexander Bernardi
 *
 * License:
 *     MIT
 */
module geo3.segment;

import geo3.point : Point3;
import geo3.scalar : isGeoScalar;


/**
 * A line segment between two points in a three-dimensional Euclidean space.
 *
 * Supported scalar types are `int`, `long`, `float`, `double`, and `real`.
 *
 * `Segment3.init` is the degenerate segment from the origin to the origin.
 *
 * Endpoint order is part of the stored value representation, but does not
 * imply traversal direction. Reversing the endpoints therefore produces a
 * different value unless both endpoints are equal.
 *
 * Degenerate segments with equal endpoints are valid.
 *
 * Equality is exact endpoint equality according to the equality semantics
 * of `Point3!T`; no tolerance or epsilon is applied.
 *
 * Floating-point endpoints may contain non-finite coordinates. `isFinite`
 * reports whether both endpoints contain only finite coordinates.
 */
struct Segment3(T)
if (isGeoScalar!T)
{
private:
    Point3!T _a;
    Point3!T _b;

public:
    /**
     * Constructs a segment from its endpoints.
     */
    this(Point3!T a, Point3!T b)
        pure nothrow @safe @nogc
    {
        _a = a;
        _b = b;
    }


    /// First stored endpoint.
    @property Point3!T a() const
        pure nothrow @safe @nogc
    {
        return _a;
    }


    /// Second stored endpoint.
    @property Point3!T b() const
        pure nothrow @safe @nogc
    {
        return _b;
    }


    /**
     * True when both endpoints contain only finite coordinates.
     */
    @property bool isFinite() const
        pure nothrow @safe @nogc
    {
        return _a.isFinite && _b.isFinite;
    }
}


/// Example constructing and inspecting a segment.
@safe unittest
{
    import geo3;

    alias P = Point3!double;
    alias S = Segment3!double;

    auto segment =
        S(
            P(1.0, 2.0, 3.0),
            P(4.0, 6.0, 8.0)
        );

    assert(
        segment.a ==
        P(1.0, 2.0, 3.0)
    );

    assert(
        segment.b ==
        P(4.0, 6.0, 8.0)
    );

    assert(segment.isFinite);
}


@safe unittest
{
    import std.meta : AliasSeq;

    /*
     * Segment3.init is the degenerate origin segment for every
     * supported scalar.
     */
    static foreach (T; AliasSeq!(int, long, float, double, real))
    {
        static assert(
            Segment3!T.init.a ==
            Point3!T.init
        );

        static assert(
            Segment3!T.init.b ==
            Point3!T.init
        );

        static assert(
            Segment3!T.init.isFinite
        );
    }


    /*
     * Unsupported scalar types must not instantiate Segment3.
     */
    static assert(!__traits(compiles, Segment3!byte));
    static assert(!__traits(compiles, Segment3!short));
    static assert(!__traits(compiles, Segment3!uint));
    static assert(!__traits(compiles, Segment3!ulong));


    alias P = Point3!double;
    alias S = Segment3!double;

    auto p = P(1.0, 2.0, 3.0);
    auto q = P(4.0, 6.0, 8.0);

    auto segment = S(p, q);

    assert(segment.a == p);
    assert(segment.b == q);
    assert(segment.isFinite);


    /*
     * Degenerate segments are valid values.
     */
    auto degenerate = S(p, p);

    assert(degenerate.a == p);
    assert(degenerate.b == p);
    assert(degenerate.isFinite);


    /*
     * Equality is exact stored-value equality.
     *
     * Reversing endpoints does not produce the same stored value.
     */
    assert(segment == S(p, q));
    assert(segment != S(q, p));


    /*
     * Non-finite endpoints remain representable.
     */
    auto infinite = S(
        P(double.infinity, 0.0, 0.0),
        P(1.0, 2.0, 3.0)
    );

    assert(!infinite.isFinite);

    auto nanSegment = S(
        P(0.0, 0.0, double.nan),
        P(1.0, 2.0, 3.0)
    );

    assert(!nanSegment.isFinite);
    assert(nanSegment != nanSegment);
}
