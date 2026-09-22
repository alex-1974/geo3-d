/**
 * Three-dimensional Euclidean point primitives.
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
module geo3.point;

import geo3.scalar : isGeoScalar;
import geo3.vector : Vector3;

import mathTraits = std.math.traits;
import std.meta : AliasSeq;
import std.traits : isFloatingPoint;


/**
 * A position in a three-dimensional Euclidean affine space.
 *
 * Supported scalar types are `int`, `long`, `float`, `double`, and `real`.
 *
 * `Point3.init` is the origin `(0, 0, 0)`.
 *
 * `Point3` follows affine point semantics. A point may be translated by
 * adding or subtracting a vector, and subtracting two points produces the
 * vector between them. Point-point addition, unary point negation, and
 * scalar multiplication or division of points are deliberately unavailable.
 *
 * Equality is exact and component-wise according to the equality semantics
 * of `T`; no tolerance or epsilon is applied. Consequently, a floating-point
 * point containing NaN does not compare equal to itself.
 *
 * Floating-point coordinates may be non-finite. `isFinite` reports whether
 * all three coordinates are finite. Integral points are always finite.
 */
struct Point3(T)
if (isGeoScalar!T)
{
private:
    // D floating-point fields default to NaN. geo3-d defines
    // Point3.init explicitly as the origin.
    T _x = 0;
    T _y = 0;
    T _z = 0;

public:
    /**
     * Constructs a point from its coordinates.
     */
    this(T x, T y, T z)
        pure nothrow @safe @nogc
    {
        _x = x;
        _y = y;
        _z = z;
    }


    /// X coordinate.
    @property T x() const
        pure nothrow @safe @nogc
    {
        return _x;
    }


    /// Y coordinate.
    @property T y() const
        pure nothrow @safe @nogc
    {
        return _y;
    }


    /// Z coordinate.
    @property T z() const
        pure nothrow @safe @nogc
    {
        return _z;
    }


    /**
     * Returns true when all three coordinates are finite.
     *
     * Integral points are always finite.
     */
    @property bool isFinite() const
        pure nothrow @safe @nogc
    {
        static if (isFloatingPoint!T)
            return mathTraits.isFinite(_x)
                && mathTraits.isFinite(_y)
                && mathTraits.isFinite(_z);
        else
            return true;
    }


    /// Translates this point by a vector.
    Point3 opBinary(string op : "+")(Vector3!T rhs)
        const pure nothrow @safe @nogc
    {
        return Point3(
            _x + rhs.x,
            _y + rhs.y,
            _z + rhs.z
        );
    }


    /// Supports Vector + Point.
    Point3 opBinaryRight(string op : "+")(Vector3!T lhs)
        const pure nothrow @safe @nogc
    {
        return Point3(
            lhs.x + _x,
            lhs.y + _y,
            lhs.z + _z
        );
    }


    /// Translates this point by the inverse of a vector.
    Point3 opBinary(string op : "-")(Vector3!T rhs)
        const pure nothrow @safe @nogc
    {
        return Point3(
            _x - rhs.x,
            _y - rhs.y,
            _z - rhs.z
        );
    }


    /// Difference between two points.
    Vector3!T opBinary(string op : "-")(Point3 rhs)
        const pure nothrow @safe @nogc
    {
        return Vector3!T(
            _x - rhs._x,
            _y - rhs._y,
            _z - rhs._z
        );
    }


    /// Compound point translation.
    void opOpAssign(string op : "+")(Vector3!T rhs)
        pure nothrow @safe @nogc
    {
        _x += rhs.x;
        _y += rhs.y;
        _z += rhs.z;
    }


    /// Compound inverse point translation.
    void opOpAssign(string op : "-")(Vector3!T rhs)
        pure nothrow @safe @nogc
    {
        _x -= rhs.x;
        _y -= rhs.y;
        _z -= rhs.z;
    }
}


@safe unittest
{
    import geo3;

    alias P = Point3!double;
    alias V = Vector3!double;

    auto point =
        P(1.0, 2.0, 3.0);

    auto shift =
        V(3.0, -1.0, 2.0);

    assert(point.x == 1.0);
    assert(point.y == 2.0);
    assert(point.z == 3.0);
    assert(point.isFinite);

    assert(
        point + shift ==
        P(4.0, 1.0, 5.0)
    );

    assert(
        shift + point ==
        P(4.0, 1.0, 5.0)
    );

    point += shift;

    assert(
        point ==
        P(4.0, 1.0, 5.0)
    );
}


@safe unittest
{
    static assert(!__traits(compiles, Point3!byte));
    static assert(!__traits(compiles, Point3!short));
    static assert(!__traits(compiles, Point3!uint));
    static assert(!__traits(compiles, Point3!ulong));

    static foreach (T; AliasSeq!(int, long, float, double, real))
    {
        static assert(Point3!T.init.x == T(0));
        static assert(Point3!T.init.y == T(0));
        static assert(Point3!T.init.z == T(0));

        static assert(
            Point3!T.init ==
            Point3!T(T(0), T(0), T(0))
        );
    }

    alias P = Point3!double;
    alias V = Vector3!double;

    auto p = P(1.0, 2.0, 3.0);
    auto q = P(4.0, 6.0, 8.0);
    auto v = V(3.0, 4.0, 5.0);

    assert(p + v == q);
    assert(v + p == q);
    assert(q - v == p);
    assert(q - p == v);

    auto r = p;

    r += v;
    assert(r == q);

    r -= v;
    assert(r == p);

    /*
     * Forbidden affine algebra is part of the API contract.
     */
    static assert(!__traits(compiles, p + q));
    static assert(!__traits(compiles, -p));
    static assert(!__traits(compiles, p * 2.0));
    static assert(!__traits(compiles, 2.0 * p));
    static assert(!__traits(compiles, p / 2.0));

    /*
     * Mixed geometry scalar types are deliberately unavailable.
     */
    Point3!int ip;
    Vector3!double dv;

    static assert(!__traits(compiles, ip + dv));

    assert(Point3!int.init.isFinite);
    assert(Point3!double.init.isFinite);

    auto nanPoint =
        P(double.nan, 0.0, 0.0);

    assert(!nanPoint.isFinite);
    assert(nanPoint != nanPoint);

    auto infinitePoint =
        P(double.infinity, 0.0, 0.0);

    assert(!infinitePoint.isFinite);

    auto nonFiniteZPoint =
        P(0.0, 0.0, double.nan);

    assert(!nonFiniteZPoint.isFinite);
}
