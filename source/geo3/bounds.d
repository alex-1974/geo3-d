/**
 * Axis-aligned three-dimensional bounds.
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
module geo3.bounds;

import geo3.point : Point3;
import geo3.scalar : isGeoScalar;

import std.traits : isFloatingPoint;


/*
 * True when a scalar is NaN.
 *
 * Integral geo3-d scalars cannot be NaN.
 */
private bool isNaNScalar(T)(T value)
    pure nothrow @safe @nogc
if (isGeoScalar!T)
{
    static if (isFloatingPoint!T)
        return value != value;
    else
        return false;
}


/*
 * True when any coordinate is NaN.
 */
private bool hasNaN(T)(Point3!T point)
    pure nothrow @safe @nogc
if (isGeoScalar!T)
{
    return isNaNScalar(point.x)
        || isNaNScalar(point.y)
        || isNaNScalar(point.z);
}


/**
 * A closed axis-aligned bounds in a three-dimensional Euclidean space.
 *
 * Supported scalar types are `int`, `long`, `float`, `double`, and `real`.
 *
 * `Bounds3.init` is empty.
 *
 * Bounds3 has two semantic states:
 *
 * - empty: contains no point;
 * - non-empty: min <= max component-wise.
 *
 * A degenerate bounds with min == max is non-empty.
 *
 * Floating-point non-empty bounds may contain infinities but never NaN.
 *
 * All public Bounds3 value operations perform no allocation and run in
 * O(1) time with O(1) auxiliary space.
 */
struct Bounds3(T)
if (isGeoScalar!T)
{
private:
    /*
     * The explicit state is intentional.
     *
     * Point3.init is the origin, while Bounds3.init must be empty.
     * Empty-state representation is not inferred from coordinate
     * sentinels such as NaN or reversed extrema.
     */
    Point3!T _min;
    Point3!T _max;
    bool _hasValue;

public:
    /**
     * True when this bounds contains no point.
     */
    @property bool empty() const
        pure nothrow @safe @nogc
    {
        return !_hasValue;
    }


    /**
     * Minimum corner.
     *
     * Precondition:
     *     This bounds is not empty.
     */
    @property Point3!T min() const
        pure nothrow @safe @nogc
    {
        assert(
            _hasValue,
            "Bounds3.min is undefined for an empty bounds"
        );

        return _min;
    }


    /**
     * Maximum corner.
     *
     * Precondition:
     *     This bounds is not empty.
     */
    @property Point3!T max() const
        pure nothrow @safe @nogc
    {
        assert(
            _hasValue,
            "Bounds3.max is undefined for an empty bounds"
        );

        return _max;
    }


    /**
     * True when all represented coordinates are finite.
     *
     * Empty bounds are finite vacuously.
     *
     * Infinite coordinates are valid Bounds3 coordinates, but cause
     * this property to return false.
     */
    @property bool isFinite() const
        pure nothrow @safe @nogc
    {
        if (!_hasValue)
            return true;

        return _min.isFinite && _max.isFinite;
    }


    /**
     * Constructs a degenerate non-empty bounds containing exactly p.
     *
     * Returns false when p contains NaN.
     *
     * On failure, result is Bounds3.init.
     */
    static bool tryFromPoint(
        Point3!T p,
        out Bounds3 result
    )
        pure nothrow @safe @nogc
    {
        if (hasNaN(p))
            return false;

        result._min = p;
        result._max = p;
        result._hasValue = true;

        return true;
    }


    /**
     * Constructs a non-empty bounds from minimum and maximum corners.
     *
     * Returns false when:
     *
     * - either corner contains NaN; or
     * - min is greater than max on any axis.
     *
     * Infinities are permitted when the ordering invariant holds.
     *
     * On failure, result is Bounds3.init.
     */
    static bool tryFromMinMax(
        Point3!T minimum,
        Point3!T maximum,
        out Bounds3 result
    )
        pure nothrow @safe @nogc
    {
        if (hasNaN(minimum) || hasNaN(maximum))
            return false;

        if (
            minimum.x > maximum.x ||
            minimum.y > maximum.y ||
            minimum.z > maximum.z
        )
        {
            return false;
        }

        result._min = minimum;
        result._max = maximum;
        result._hasValue = true;

        return true;
    }


    /**
     * Extends this bounds so that it contains p.
     *
     * Returns false when p contains NaN.
     *
     * On failure, this bounds remains unchanged.
     */
    bool tryExtend(Point3!T p)
        pure nothrow @safe @nogc
    {
        if (hasNaN(p))
            return false;

        if (!_hasValue)
        {
            _min = p;
            _max = p;
            _hasValue = true;

            return true;
        }

        const T minX =
            p.x < _min.x ? p.x : _min.x;

        const T minY =
            p.y < _min.y ? p.y : _min.y;

        const T minZ =
            p.z < _min.z ? p.z : _min.z;

        const T maxX =
            p.x > _max.x ? p.x : _max.x;

        const T maxY =
            p.y > _max.y ? p.y : _max.y;

        const T maxZ =
            p.z > _max.z ? p.z : _max.z;

        _min = Point3!T(
            minX,
            minY,
            minZ
        );

        _max = Point3!T(
            maxX,
            maxY,
            maxZ
        );

        return true;
    }


    /**
     * Extends this bounds so that it contains other.
     *
     * Extending by an empty bounds has no effect.
     */
    void extend(Bounds3 other)
        pure nothrow @safe @nogc
    {
        if (!other._hasValue)
            return;

        if (!_hasValue)
        {
            _min = other._min;
            _max = other._max;
            _hasValue = true;

            return;
        }

        const T minX =
            other._min.x < _min.x
                ? other._min.x
                : _min.x;

        const T minY =
            other._min.y < _min.y
                ? other._min.y
                : _min.y;

        const T minZ =
            other._min.z < _min.z
                ? other._min.z
                : _min.z;

        const T maxX =
            other._max.x > _max.x
                ? other._max.x
                : _max.x;

        const T maxY =
            other._max.y > _max.y
                ? other._max.y
                : _max.y;

        const T maxZ =
            other._max.z > _max.z
                ? other._max.z
                : _max.z;

        _min = Point3!T(
            minX,
            minY,
            minZ
        );

        _max = Point3!T(
            maxX,
            maxY,
            maxZ
        );
    }


    /**
     * True when this closed bounds contains p.
     *
     * Empty bounds contain no point.
     * A point containing NaN is never contained.
     */
    bool contains(Point3!T p) const
        pure nothrow @safe @nogc
    {
        if (!_hasValue || hasNaN(p))
            return false;

        return p.x >= _min.x
            && p.x <= _max.x
            && p.y >= _min.y
            && p.y <= _max.y
            && p.z >= _min.z
            && p.z <= _max.z;
    }


    /**
     * True when this closed bounds intersects other.
     *
     * Empty bounds never intersect.
     *
     * Touching faces, edges, or corners count as intersection.
     */
    bool intersects(Bounds3 other) const
        pure nothrow @safe @nogc
    {
        if (!_hasValue || !other._hasValue)
            return false;

        return !(
               _max.x < other._min.x
            || other._max.x < _min.x
            || _max.y < other._min.y
            || other._max.y < _min.y
            || _max.z < other._min.z
            || other._max.z < _min.z
        );
    }


    /**
     * Exact bounds equality.
     *
     * All empty Bounds3 values of the same type compare equal
     * independently of their internal representation.
     */
    bool opEquals(const Bounds3 rhs) const
        pure nothrow @safe @nogc
    {
        if (_hasValue != rhs._hasValue)
            return false;

        if (!_hasValue)
            return true;

        return _min == rhs._min
            && _max == rhs._max;
    }
}


/// Example using the core Bounds3 value operations.
@safe unittest
{
    import geo3;

    alias P = Point3!int;
    alias B = Bounds3!int;

    assert(B.init.empty);
    assert(B.init.isFinite);

    B first;

    assert(
        B.tryFromPoint(
            P(3, 4, 5),
            first
        )
    );

    assert(!first.empty);
    assert(first.min == P(3, 4, 5));
    assert(first.max == P(3, 4, 5));
    assert(first.isFinite);

    B second;

    assert(
        B.tryFromPoint(
            P(-2, 8, 1),
            second
        )
    );

    first.extend(second);

    assert(first.min == P(-2, 4, 1));
    assert(first.max == P(3, 8, 5));

    const copy = first;
    assert(copy == first);
}


@safe unittest
{
    import std.meta : AliasSeq;

    /*
     * Bounds3.init is empty for every supported scalar.
     */
    static foreach (T; AliasSeq!(int, long, float, double, real))
    {
        static assert(Bounds3!T.init.empty);
        static assert(Bounds3!T.init.isFinite);
    }


    /*
     * Unsupported scalar types must remain unavailable.
     */
    static assert(!__traits(compiles, Bounds3!byte));
    static assert(!__traits(compiles, Bounds3!short));
    static assert(!__traits(compiles, Bounds3!uint));
    static assert(!__traits(compiles, Bounds3!ulong));


    alias P = Point3!double;
    alias B = Bounds3!double;


    /*
     * Empty and degenerate bounds are distinct.
     */
    B pointBounds;

    assert(
        B.tryFromPoint(
            P(2.0, 3.0, 4.0),
            pointBounds
        )
    );

    assert(!pointBounds.empty);
    assert(pointBounds.min == P(2.0, 3.0, 4.0));
    assert(pointBounds.max == P(2.0, 3.0, 4.0));
    assert(pointBounds != B.init);


    /*
     * Ordered min/max construction.
     */
    B bounds;

    assert(
        B.tryFromMinMax(
            P(1.0, 2.0, 3.0),
            P(4.0, 6.0, 8.0),
            bounds
        )
    );

    assert(!bounds.empty);
    assert(bounds.min == P(1.0, 2.0, 3.0));
    assert(bounds.max == P(4.0, 6.0, 8.0));


    /*
     * Reversal on every individual axis is rejected.
     */
    B invalid;

    assert(
        !B.tryFromMinMax(
            P(5.0, 2.0, 3.0),
            P(4.0, 6.0, 8.0),
            invalid
        )
    );

    assert(invalid.empty);

    assert(
        !B.tryFromMinMax(
            P(1.0, 7.0, 3.0),
            P(4.0, 6.0, 8.0),
            invalid
        )
    );

    assert(invalid.empty);

    assert(
        !B.tryFromMinMax(
            P(1.0, 2.0, 9.0),
            P(4.0, 6.0, 8.0),
            invalid
        )
    );

    assert(invalid.empty);


    /*
     * NaN is rejected on every coordinate family and must never become
     * an empty-state synonym.
     */
    B nanBounds;

    assert(
        !B.tryFromPoint(
            P(0.0, 0.0, double.nan),
            nanBounds
        )
    );

    assert(nanBounds.empty);

    assert(
        !B.tryFromMinMax(
            P(0.0, double.nan, 0.0),
            P(1.0, 1.0, 1.0),
            nanBounds
        )
    );

    assert(nanBounds.empty);


    /*
     * Infinity is valid when ordering remains meaningful.
     */
    B infiniteBounds;

    assert(
        B.tryFromMinMax(
            P(
                -double.infinity,
                0.0,
                -double.infinity
            ),
            P(
                double.infinity,
                10.0,
                double.infinity
            ),
            infiniteBounds
        )
    );

    assert(!infiniteBounds.empty);
    assert(!infiniteBounds.isFinite);

    assert(
        infiniteBounds.contains(
            P(0.0, 5.0, 100.0)
        )
    );

    assert(
        infiniteBounds.contains(
            P(
                double.infinity,
                5.0,
                double.infinity
            )
        )
    );


    /*
     * Incremental extension from empty.
     */
    B accumulated;

    assert(
        accumulated.tryExtend(
            P(5.0, 4.0, 9.0)
        )
    );

    assert(
        accumulated.tryExtend(
            P(1.0, 7.0, 2.0)
        )
    );

    assert(
        accumulated.min ==
        P(1.0, 4.0, 2.0)
    );

    assert(
        accumulated.max ==
        P(5.0, 7.0, 9.0)
    );


    /*
     * Failed point extension preserves the original value.
     */
    const before = accumulated;

    assert(
        !accumulated.tryExtend(
            P(10.0, 10.0, double.nan)
        )
    );

    assert(accumulated == before);


    /*
     * Bounds extension.
     */
    B other;

    assert(
        B.tryFromMinMax(
            P(-2.0, 5.0, -1.0),
            P(3.0, 9.0, 12.0),
            other
        )
    );

    accumulated.extend(other);

    assert(
        accumulated.min ==
        P(-2.0, 4.0, -1.0)
    );

    assert(
        accumulated.max ==
        P(5.0, 9.0, 12.0)
    );

    auto unchanged = accumulated;
    unchanged.extend(B.init);

    assert(unchanged == accumulated);


    /*
     * Closed containment semantics, including Z.
     */
    assert(
        !B.init.contains(
            P(0.0, 0.0, 0.0)
        )
    );

    assert(bounds.contains(P(1.0, 2.0, 3.0)));
    assert(bounds.contains(P(4.0, 6.0, 8.0)));
    assert(bounds.contains(P(2.0, 3.0, 5.0)));

    assert(!bounds.contains(P(2.0, 3.0, 9.0)));
    assert(!bounds.contains(P(2.0, 3.0, double.nan)));


    /*
     * Closed intersection semantics.
     *
     * Touching on a Z face counts as intersection.
     */
    B touching;
    B separate;

    assert(
        B.tryFromMinMax(
            P(2.0, 3.0, 8.0),
            P(3.0, 5.0, 10.0),
            touching
        )
    );

    assert(
        B.tryFromMinMax(
            P(2.0, 3.0, 8.0001),
            P(3.0, 5.0, 10.0),
            separate
        )
    );

    assert(bounds.intersects(touching));
    assert(touching.intersects(bounds));

    assert(!bounds.intersects(separate));
    assert(!separate.intersects(bounds));

    assert(!bounds.intersects(B.init));
    assert(!B.init.intersects(bounds));


    /*
     * Equality.
     */
    assert(B.init == B.init);
    assert(bounds == bounds);
    assert(bounds != pointBounds);
}
