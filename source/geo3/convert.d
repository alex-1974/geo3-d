/**
 * Checked scalar conversion of three-dimensional geometry values.
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
module geo3.convert;

import geo3.point :
    Point3;

import geo3.scalar :
    isGeoScalar;

import geo3.segment :
    Segment3;

import geo3.vector :
    Vector3;

import std.math.rounding :
    ceil,
    floor;

import std.math.traits :
    isFinite;

import std.traits :
    isFloatingPoint;


private enum bool isGeoIntegral(T) =
       is(T == int)
    || is(T == long);


/*
 * Removes the fractional part while retaining the floating-point type.
 */
private T truncateTowardZero(T)(T value)
    pure nothrow @safe @nogc
if (isFloatingPoint!T)
{
    if (value > T(0))
        return floor(value);

    if (value < T(0))
        return ceil(value);

    /*
     * Preserves zero, signed zero and NaN.
     */
    return value;
}


/*
 * Checked scalar conversion used by the public geometry conversions.
 *
 * Precision loss inside the representable target range is permitted.
 * Range overflow is not.
 *
 * Floating-point to integral conversion additionally requires:
 *
 * - a finite source;
 * - an already integral-valued source;
 * - a value inside the target range.
 */
private bool tryScalarConvert(To, From)(
    From value,
    out To result
)
    pure nothrow @safe @nogc
if (isGeoScalar!To && isGeoScalar!From)
{
    static if (
        isGeoIntegral!From &&
        isGeoIntegral!To
    )
    {
        static if (From.sizeof <= To.sizeof)
        {
            result = cast(To) value;
            return true;
        }
        else
        {
            if (
                value < cast(From) To.min ||
                value > cast(From) To.max
            )
            {
                return false;
            }

            result = cast(To) value;
            return true;
        }
    }
    else static if (
        isGeoIntegral!From &&
        isFloatingPoint!To
    )
    {
        /*
         * int and long fit inside the exponent range of every supported
         * floating-point type. Precision loss is permitted.
         */
        result = cast(To) value;
        return true;
    }
    else static if (
        isFloatingPoint!From &&
        isGeoIntegral!To
    )
    {
        if (!isFinite(value))
            return false;

        if (value != truncateTowardZero(value))
            return false;

        /*
         * Avoid comparing against cast(From) To.max.
         *
         * For example, cast(double) long.max rounds to 2^63, which is
         * already outside the positive long range.
         *
         * For a signed N-bit integer the valid interval is:
         *
         *     [-2^(N-1), 2^(N-1))
         *
         * Powers of two are exactly representable in binary floating point.
         */
        enum shift =
            To.sizeof * 8 - 1;

        enum ulong upperMagnitude =
            1UL << shift;

        const From upperExclusive =
            cast(From) upperMagnitude;

        const From lowerInclusive =
            -upperExclusive;

        if (
            value < lowerInclusive ||
            value >= upperExclusive
        )
        {
            return false;
        }

        result = cast(To) value;
        return true;
    }
    else static if (
        isFloatingPoint!From &&
        isFloatingPoint!To
    )
    {
        const To converted =
            cast(To) value;

        /*
         * NaN and ±infinity remain valid floating-point coordinates.
         *
         * A finite source becoming infinity in the target type is range
         * overflow and therefore rejected.
         */
        if (
            isFinite(value) &&
            !isFinite(converted)
        )
        {
            return false;
        }

        result = converted;
        return true;
    }
    else
    {
        static assert(
            0,
            "unsupported geo3-d scalar conversion"
        );
    }
}


/**
 * Checked component-wise conversion of a point.
 *
 * Source and target scalar types must belong to the geo3-d scalar domain.
 *
 * Conversion is explicit and checked for range, but precision loss within
 * the representable target range is permitted.
 *
 * In particular:
 *
 * - integral-to-integral narrowing fails outside the target range;
 * - integral-to-floating conversion may lose precision;
 * - floating-to-integral conversion requires every coordinate to be finite,
 *   already integral-valued, and inside the target integral range;
 * - floating-to-floating conversion permits NaN and infinity;
 * - floating-to-floating conversion fails when a finite source would become
 *   infinity in the target type.
 *
 * No implicit rounding or quantisation is performed.
 *
 * Returns false when any coordinate cannot be converted according to these
 * rules.
 *
 * On failure, `result` is `Point3!To.init`.
 *
 * No allocation is performed.
 *
 * Complexity:
 *     O(1) time and O(1) auxiliary space.
 */
bool tryConvert(To, From)(
    Point3!From source,
    out Point3!To result
)
    pure nothrow @safe @nogc
if (isGeoScalar!To && isGeoScalar!From)
{
    To x;
    To y;
    To z;

    if (!tryScalarConvert!To(source.x, x))
        return false;

    if (!tryScalarConvert!To(source.y, y))
        return false;

    if (!tryScalarConvert!To(source.z, z))
        return false;

    result = Point3!To(
        x,
        y,
        z
    );

    return true;
}


/**
 * Checked component-wise conversion of a vector.
 *
 * Conversion follows exactly the same scalar rules as Point3 conversion:
 * range overflow is rejected, precision loss inside the representable
 * target range is permitted, and no implicit rounding is performed.
 *
 * Returns false when any component cannot be converted.
 *
 * On failure, `result` is `Vector3!To.init`.
 *
 * No allocation is performed.
 *
 * Complexity:
 *     O(1) time and O(1) auxiliary space.
 */
bool tryConvert(To, From)(
    Vector3!From source,
    out Vector3!To result
)
    pure nothrow @safe @nogc
if (isGeoScalar!To && isGeoScalar!From)
{
    To x;
    To y;
    To z;

    if (!tryScalarConvert!To(source.x, x))
        return false;

    if (!tryScalarConvert!To(source.y, y))
        return false;

    if (!tryScalarConvert!To(source.z, z))
        return false;

    result = Vector3!To(
        x,
        y,
        z
    );

    return true;
}


/**
 * Checked component-wise conversion of a segment.
 *
 * Both endpoints are converted using the Point3 conversion rules. The
 * complete operation succeeds only when both endpoints convert successfully.
 *
 * Precision loss inside the representable target range is permitted.
 * No implicit rounding or quantisation is performed.
 *
 * On failure, `result` is `Segment3!To.init`.
 *
 * No allocation is performed.
 *
 * Complexity:
 *     O(1) time and O(1) auxiliary space.
 */
bool tryConvert(To, From)(
    Segment3!From source,
    out Segment3!To result
)
    pure nothrow @safe @nogc
if (isGeoScalar!To && isGeoScalar!From)
{
    Point3!To a;
    Point3!To b;

    if (!tryConvert!To(source.a, a))
        return false;

    if (!tryConvert!To(source.b, b))
        return false;

    result = Segment3!To(
        a,
        b
    );

    return true;
}


/// Example of checked three-dimensional conversion without implicit rounding.
@safe unittest
{
    import geo3;

    Point3!int result;

    assert(
        Point3!double(
            3.0,
            -2.0,
            7.0
        ).tryConvert(result)
    );

    assert(
        result ==
        Point3!int(
            3,
            -2,
            7
        )
    );

    /*
     * Failure on Z must invalidate the complete conversion.
     */
    assert(
        !Point3!double(
            3.0,
            -2.0,
            7.5
        ).tryConvert(result)
    );

    assert(
        result ==
        Point3!int.init
    );
}


// Regression coverage.
@safe unittest
{
    /*
     * Integral widening.
     */
    Point3!long pointLong;

    assert(
        Point3!int(
            1,
            -2,
            3
        ).tryConvert(pointLong)
    );

    assert(
        pointLong ==
        Point3!long(
            1,
            -2,
            3
        )
    );


    /*
     * Integral narrowing is checked independently on every axis.
     */
    Point3!int pointInt;

    assert(
        Point3!long(
            123,
            -456,
            789
        ).tryConvert(pointInt)
    );

    assert(
        pointInt ==
        Point3!int(
            123,
            -456,
            789
        )
    );

    assert(
        !Point3!long(
            0,
            0,
            long.max
        ).tryConvert(pointInt)
    );

    assert(
        pointInt ==
        Point3!int.init
    );


    /*
     * Integral to floating conversion is range-valid. Precision loss is
     * explicitly permitted.
     */
    Point3!double pointDouble;

    assert(
        Point3!long(
            123,
            -456,
            long.max
        ).tryConvert(pointDouble)
    );

    assert(pointDouble.x == 123.0);
    assert(pointDouble.y == -456.0);
    assert(pointDouble.z > 0.0);


    /*
     * Floating to integral requires already integral-valued coordinates.
     */
    assert(
        Point3!double(
            3.0,
            -4.0,
            5.0
        ).tryConvert(pointInt)
    );

    assert(
        pointInt ==
        Point3!int(
            3,
            -4,
            5
        )
    );

    assert(
        !Point3!double(
            3.0,
            -4.0,
            5.25
        ).tryConvert(pointInt)
    );

    assert(
        pointInt ==
        Point3!int.init
    );


    /*
     * Non-finite floating-point values cannot become integers.
     * Exercise all three coordinate positions.
     */
    assert(
        !Point3!double(
            double.nan,
            0.0,
            0.0
        ).tryConvert(pointInt)
    );

    assert(
        pointInt ==
        Point3!int.init
    );

    assert(
        !Point3!double(
            0.0,
            double.infinity,
            0.0
        ).tryConvert(pointInt)
    );

    assert(
        pointInt ==
        Point3!int.init
    );

    assert(
        !Point3!double(
            0.0,
            0.0,
            -double.infinity
        ).tryConvert(pointInt)
    );

    assert(
        pointInt ==
        Point3!int.init
    );


    /*
     * Floating-point non-finite values remain representable when converting
     * to another floating-point geometry type.
     */
    Point3!float pointFloat;

    assert(
        Point3!double(
            double.infinity,
            1.0,
            -double.infinity
        ).tryConvert(pointFloat)
    );

    assert(
        pointFloat.x ==
        float.infinity
    );

    assert(
        pointFloat.z ==
        -float.infinity
    );

    assert(
        Point3!double(
            1.0,
            double.nan,
            2.0
        ).tryConvert(pointFloat)
    );

    assert(
        pointFloat.y !=
        pointFloat.y
    );


    /*
     * Finite floating-point narrowing must not overflow to infinity.
     * Failure on Z also confirms that the third coordinate participates.
     */
    assert(
        !Point3!double(
            0.0,
            0.0,
            double.max
        ).tryConvert(pointFloat)
    );

    assert(
        pointFloat ==
        Point3!float.init
    );


    /*
     * The upper boundary of long is intentionally tested using 2^63.
     *
     * cast(double) long.max is also 2^63, so comparing a double directly
     * against cast(double) long.max would incorrectly accept this value.
     */
    Point3!long longBoundary;

    assert(
        !Point3!double(
            0.0,
            0.0,
            9_223_372_036_854_775_808.0
        ).tryConvert(longBoundary)
    );

    assert(
        longBoundary ==
        Point3!long.init
    );

    /*
     * -2^63 is valid and exactly representable.
     */
    assert(
        Point3!double(
            0.0,
            0.0,
            -9_223_372_036_854_775_808.0
        ).tryConvert(longBoundary)
    );

    assert(
        longBoundary.z ==
        long.min
    );


    /*
     * Vector conversion follows the same scalar rules.
     */
    Vector3!int vectorInt;

    assert(
        Vector3!double(
            4.0,
            -7.0,
            12.0
        ).tryConvert(vectorInt)
    );

    assert(
        vectorInt ==
        Vector3!int(
            4,
            -7,
            12
        )
    );

    assert(
        !Vector3!double(
            4.0,
            -7.0,
            12.25
        ).tryConvert(vectorInt)
    );

    assert(
        vectorInt ==
        Vector3!int.init
    );


    /*
     * Segment conversion reuses the Point3 rules for both endpoints.
     */
    Segment3!int segmentInt;

    assert(
        Segment3!double(
            Point3!double(
                1.0,
                -2.0,
                3.0
            ),
            Point3!double(
                4.0,
                5.0,
                -6.0
            )
        ).tryConvert(segmentInt)
    );

    assert(
        segmentInt ==
        Segment3!int(
            Point3!int(
                1,
                -2,
                3
            ),
            Point3!int(
                4,
                5,
                -6
            )
        )
    );


    /*
     * Failure in the first endpoint is transactional.
     */
    assert(
        !Segment3!double(
            Point3!double(
                1.0,
                -2.0,
                3.5
            ),
            Point3!double(
                4.0,
                5.0,
                6.0
            )
        ).tryConvert(segmentInt)
    );

    assert(
        segmentInt ==
        Segment3!int.init
    );


    /*
     * Failure in the second endpoint, including its Z component, is equally
     * transactional.
     */
    assert(
        !Segment3!double(
            Point3!double(
                1.0,
                -2.0,
                3.0
            ),
            Point3!double(
                4.0,
                5.0,
                6.5
            )
        ).tryConvert(segmentInt)
    );

    assert(
        segmentInt ==
        Segment3!int.init
    );
}
