/**
 * Three-dimensional Euclidean vector primitives.
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
module geo3.vector;

import geo3.scalar : isGeoScalar;

import mathTraits = std.math.traits;
import std.meta : AliasSeq;
import std.traits : isFloatingPoint;


/**
 * A displacement in a three-dimensional Euclidean vector space.
 *
 * Supported scalar types are `int`, `long`, `float`, `double`, and `real`.
 *
 * `Vector3.init` is the zero vector `(0, 0, 0)`.
 *
 * Vectors support addition, subtraction, negation, and scalar
 * multiplication. Scalar division is available when the resulting scalar
 * type is floating-point. Arithmetic result types follow D's scalar
 * arithmetic rules subject to the supported `geo3-d` scalar types.
 *
 * Equality is exact and component-wise according to the equality semantics
 * of `T`; no tolerance or epsilon is applied. Consequently, a floating-point
 * vector containing NaN does not compare equal to itself.
 *
 * Floating-point components may be non-finite. `isFinite` reports whether
 * all three components are finite. Integral vectors are always finite.
 */
struct Vector3(T)
if (isGeoScalar!T)
{
private:
    // D floating-point fields default to NaN. geo3-d defines
    // Vector3.init explicitly as the zero vector.
    T _x = 0;
    T _y = 0;
    T _z = 0;

public:
    /**
     * Constructs a vector from its components.
     */
    this(T x, T y, T z)
        pure nothrow @safe @nogc
    {
        _x = x;
        _y = y;
        _z = z;
    }


    /// X component.
    @property T x() const
        pure nothrow @safe @nogc
    {
        return _x;
    }


    /// Y component.
    @property T y() const
        pure nothrow @safe @nogc
    {
        return _y;
    }


    /// Z component.
    @property T z() const
        pure nothrow @safe @nogc
    {
        return _z;
    }


    /**
     * Returns true when all three components are finite.
     *
     * Integral vectors are always finite.
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


    /// Vector addition.
    Vector3 opBinary(string op : "+")(Vector3 rhs)
        const pure nothrow @safe @nogc
    {
        return Vector3(
            _x + rhs._x,
            _y + rhs._y,
            _z + rhs._z
        );
    }


    /// Vector subtraction.
    Vector3 opBinary(string op : "-")(Vector3 rhs)
        const pure nothrow @safe @nogc
    {
        return Vector3(
            _x - rhs._x,
            _y - rhs._y,
            _z - rhs._z
        );
    }


    /// Vector negation.
    Vector3 opUnary(string op : "-")()
        const pure nothrow @safe @nogc
    {
        return Vector3(
            -_x,
            -_y,
            -_z
        );
    }


    /**
     * Scales this vector.
     *
     * The arithmetic result scalar becomes the result vector scalar.
     */
    auto opBinary(string op : "*", S)(S rhs)
        const pure nothrow @safe @nogc
    if (
        isGeoScalar!S
        && isGeoScalar!(typeof(T.init * S.init))
    )
    {
        alias R = typeof(T.init * S.init);

        return Vector3!R(
            _x * rhs,
            _y * rhs,
            _z * rhs
        );
    }


    /**
     * Scales this vector with the scalar on the left.
     */
    auto opBinaryRight(string op : "*", S)(S lhs)
        const pure nothrow @safe @nogc
    if (
        isGeoScalar!S
        && isGeoScalar!(typeof(S.init * T.init))
    )
    {
        alias R = typeof(S.init * T.init);

        return Vector3!R(
            lhs * _x,
            lhs * _y,
            lhs * _z
        );
    }


    /**
     * Divides this vector by a scalar.
     *
     * Integer-result division is deliberately unavailable.
     */
    auto opBinary(string op : "/", S)(S rhs)
        const pure nothrow @safe @nogc
    if (
        isGeoScalar!S
        && isFloatingPoint!(typeof(T.init / S.init))
    )
    {
        alias R = typeof(T.init / S.init);

        return Vector3!R(
            _x / rhs,
            _y / rhs,
            _z / rhs
        );
    }


    /// Compound vector addition.
    void opOpAssign(string op : "+")(Vector3 rhs)
        pure nothrow @safe @nogc
    {
        _x += rhs._x;
        _y += rhs._y;
        _z += rhs._z;
    }


    /// Compound vector subtraction.
    void opOpAssign(string op : "-")(Vector3 rhs)
        pure nothrow @safe @nogc
    {
        _x -= rhs._x;
        _y -= rhs._y;
        _z -= rhs._z;
    }


    /**
     * Compound scalar multiplication.
     *
     * Available only when the result scalar remains T.
     */
    void opOpAssign(string op : "*", S)(S rhs)
        pure nothrow @safe @nogc
    if (
        isGeoScalar!S
        && is(typeof(T.init * S.init) == T)
    )
    {
        _x *= rhs;
        _y *= rhs;
        _z *= rhs;
    }


    /**
     * Compound scalar division.
     *
     * Available only for floating-point vectors and when the result
     * scalar remains T.
     */
    void opOpAssign(string op : "/", S)(S rhs)
        pure nothrow @safe @nogc
    if (
        isGeoScalar!S
        && isFloatingPoint!T
        && is(typeof(T.init / S.init) == T)
    )
    {
        _x /= rhs;
        _y /= rhs;
        _z /= rhs;
    }
}


@safe unittest
{
    import geo3;

    alias V = Vector3!double;

    auto vector =
        V(2.0, -4.0, 6.0);

    assert(vector.x == 2.0);
    assert(vector.y == -4.0);
    assert(vector.z == 6.0);
    assert(vector.isFinite);

    assert(
        -vector ==
        V(-2.0, 4.0, -6.0)
    );

    assert(
        vector + V(1.0, 1.0, 1.0) ==
        V(3.0, -3.0, 7.0)
    );

    assert(
        vector * 0.5 ==
        V(1.0, -2.0, 3.0)
    );

    assert(
        0.5 * vector ==
        V(1.0, -2.0, 3.0)
    );

    vector *= 2.0;

    assert(
        vector ==
        V(4.0, -8.0, 12.0)
    );
}


@safe unittest
{
    static assert(!__traits(compiles, Vector3!byte));
    static assert(!__traits(compiles, Vector3!short));
    static assert(!__traits(compiles, Vector3!uint));
    static assert(!__traits(compiles, Vector3!ulong));

    static foreach (T; AliasSeq!(int, long, float, double, real))
    {
        static assert(Vector3!T.init.x == T(0));
        static assert(Vector3!T.init.y == T(0));
        static assert(Vector3!T.init.z == T(0));

        static assert(
            Vector3!T.init ==
            Vector3!T(T(0), T(0), T(0))
        );
    }

    alias V = Vector3!double;

    auto a = V(1.0, 2.0, 3.0);
    auto b = V(4.0, 6.0, 8.0);

    assert(a + b == V(5.0, 8.0, 11.0));
    assert(b - a == V(3.0, 4.0, 5.0));
    assert(-a == V(-1.0, -2.0, -3.0));

    assert(a * 2.0 == V(2.0, 4.0, 6.0));
    assert(2.0 * a == V(2.0, 4.0, 6.0));
    assert(a / 2.0 == V(0.5, 1.0, 1.5));

    auto c = a;

    c += b;
    assert(c == V(5.0, 8.0, 11.0));

    c -= b;
    assert(c == a);

    c *= 2;
    assert(c == V(2.0, 4.0, 6.0));

    c /= 2;
    assert(c == a);

    static assert(
        is(
            typeof(Vector3!int(2, 4, 6) * 0.5)
            == Vector3!double
        )
    );

    static assert(
        is(
            typeof(Vector3!double(2, 4, 6) * 2)
            == Vector3!double
        )
    );

    Vector3!int iv;

    static assert(!__traits(compiles, iv / 2));
    static assert(__traits(compiles, iv / 2.0));

    assert(Vector3!int.init.isFinite);
    assert(Vector3!double.init.isFinite);

    auto nanVector =
        V(double.nan, 0.0, 0.0);

    assert(!nanVector.isFinite);
    assert(nanVector != nanVector);

    auto infiniteVector =
        V(double.infinity, 0.0, 0.0);

    assert(!infiniteVector.isFinite);

    auto nonFiniteZVector =
        V(0.0, 0.0, double.nan);

    assert(!nonFiniteZVector.isFinite);
}
