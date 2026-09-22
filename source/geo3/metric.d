/**
 * Three-dimensional Euclidean metric operations.
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
module geo3.metric;

import geo3.point :
    Point3;

import geo3.polyline_view :
    Polyline3View;

import geo3.segment :
    Segment3;

import geo3.scalar :
    MetricScalar,
    isGeoScalar;

import std.math.algebraic :
    hypot;

import std.math.traits :
    isFinite;


/*
 * Supported signed integral geometry scalar types.
 */
private enum bool isMetricIntegral(T) =
       is(T == int)
    || is(T == long);


/*
 * Unsigned type capable of holding the complete magnitude range of T.
 */
private template UnsignedMetricIntegral(T)
if (isMetricIntegral!T)
{
    static if (is(T == int))
        alias UnsignedMetricIntegral = uint;
    else
        alias UnsignedMetricIntegral = ulong;
}


/*
 * Exact unsigned magnitude of a signed integral value.
 *
 * The -(value + 1) formulation avoids overflow for T.min.
 */
private UnsignedMetricIntegral!T unsignedMagnitude(T)(T value)
    pure nothrow @safe @nogc
if (isMetricIntegral!T)
{
    alias U = UnsignedMetricIntegral!T;

    if (value >= 0)
        return cast(U) value;

    return cast(U)(-(value + 1)) + U(1);
}


/*
 * Exact absolute difference between two supported signed integer values.
 *
 * No signed subtraction overflow occurs. The result may span the complete
 * corresponding unsigned type:
 *
 *     int  -> uint
 *     long -> ulong
 */
private UnsignedMetricIntegral!T unsignedDifference(T)(T a, T b)
    pure nothrow @safe @nogc
if (isMetricIntegral!T)
{
    alias U = UnsignedMetricIntegral!T;

    if ((a < 0) != (b < 0))
        return unsignedMagnitude(a) + unsignedMagnitude(b);

    if (a >= b)
        return cast(U)(a - b);

    return cast(U)(b - a);
}


/*
 * Signed component difference in the metric computation type.
 *
 * Semantics:
 *
 *     a - b
 *
 * For integral geometry the magnitude is obtained exactly before conversion
 * to MetricScalar. The sign is applied only after conversion.
 *
 * This avoids both signed integer subtraction overflow and loss of small
 * differences caused by converting large integer coordinates to floating
 * point before subtraction.
 */
private MetricScalar!T signedMetricDifference(T)(T a, T b)
    pure nothrow @safe @nogc
if (isGeoScalar!T)
{
    alias M = MetricScalar!T;

    static if (isMetricIntegral!T)
    {
        if (a == b)
            return M(0);

        const M magnitude =
            cast(M) unsignedDifference(a, b);

        return a > b
            ? magnitude
            : -magnitude;
    }
    else
    {
        return
            cast(M) a -
            cast(M) b;
    }
}


/**
 * Squared Euclidean distance between two three-dimensional points.
 *
 * Returns:
 *     The squared distance in `MetricScalar!T`.
 *
 * Integer coordinate geometry is converted to floating-point metric
 * arithmetic after each component difference has been obtained without
 * signed overflow.
 *
 * This operation does not promise exact integral arithmetic. In
 * particular, long-coordinate results may lose precision after conversion
 * to double.
 *
 * Floating-point NaN and infinity are not rejected. Results follow normal
 * floating-point arithmetic. Very large finite results may overflow to
 * infinity.
 *
 * This function is a metric computation, not a robust exact distance
 * comparison predicate.
 *
 * No allocation is performed.
 *
 * Complexity:
 *     O(1) time and O(1) auxiliary space.
 */
MetricScalar!T squaredDistance(T)(
    Point3!T a,
    Point3!T b
)
    pure nothrow @safe @nogc
if (isGeoScalar!T)
{
    alias M = MetricScalar!T;

    const M dx =
        signedMetricDifference(
            a.x,
            b.x
        );

    const M dy =
        signedMetricDifference(
            a.y,
            b.y
        );

    const M dz =
        signedMetricDifference(
            a.z,
            b.z
        );

    return
        dx * dx +
        dy * dy +
        dz * dz;
}


/// Example computing squared distance without taking a square root.
@safe unittest
{
    import geo3;

    const a =
        Point3!int(
            0,
            0,
            0
        );

    const b =
        Point3!int(
            2,
            3,
            6
        );

    static assert(
        is(
            typeof(
                squaredDistance(
                    a,
                    b
                )
            ) ==
            double
        )
    );

    assert(
        squaredDistance(
            a,
            b
        ) ==
        49.0
    );
}


/**
 * Euclidean distance between two three-dimensional points.
 *
 * Integral coordinate differences are obtained without signed overflow
 * before conversion to `MetricScalar!T`. This also preserves small
 * differences between large integer coordinates that would otherwise be
 * lost by converting the coordinates before subtraction.
 *
 * The resulting metric value is floating-point. In particular,
 * long-coordinate results may lose precision after conversion to double.
 *
 * Distance is calculated with `hypot` so avoidable intermediate square
 * overflow is not introduced.
 *
 * Floating-point NaN and infinity are not rejected. Results follow normal
 * floating-point arithmetic.
 *
 * This function is a metric computation, not a robust exact distance
 * comparison predicate.
 *
 * No allocation is performed.
 *
 * Complexity:
 *     O(1) time and O(1) auxiliary space.
 */
MetricScalar!T distance(T)(
    Point3!T a,
    Point3!T b
)
    pure nothrow @safe @nogc
if (isGeoScalar!T)
{
    const dx =
        signedMetricDifference(
            a.x,
            b.x
        );

    const dy =
        signedMetricDifference(
            a.y,
            b.y
        );

    const dz =
        signedMetricDifference(
            a.z,
            b.z
        );

    return hypot(
        hypot(
            dx,
            dy
        ),
        dz
    );
}


/// Example computing a three-dimensional Euclidean distance.
@safe unittest
{
    import geo3;

    alias P = Point3!double;

    assert(
        distance(
            P(0.0, 0.0, 0.0),
            P(2.0, 3.0, 6.0)
        ) ==
        7.0
    );
}


/**
 * Euclidean length of a segment.
 *
 * Uses the same metric computation policy as point-to-point distance.
 *
 * In particular:
 *
 * - integer coordinate differences are obtained without signed overflow;
 * - int, long and float geometry compute in double;
 * - real geometry computes in real;
 * - hypot is used indirectly through distance().
 *
 * Floating-point non-finite coordinates follow the same arithmetic
 * semantics as distance().
 *
 * This is a metric computation, not an exact topological predicate.
 *
 * No allocation is performed.
 *
 * Complexity:
 *     O(1) time and O(1) auxiliary space.
 */
MetricScalar!T segmentLength(T)(Segment3!T segment)
    pure nothrow @safe @nogc
if (isGeoScalar!T)
{
    return distance(
        segment.a,
        segment.b
    );
}


/// Example computing the Euclidean length of a 3D segment.
@safe unittest
{
    import geo3;

    const segment =
        Segment3!double(
            Point3!double(
                0.0,
                0.0,
                0.0
            ),
            Point3!double(
                2.0,
                3.0,
                6.0
            )
        );

    assert(
        segmentLength(segment) ==
        7.0
    );
}


/**
 * Euclidean length of a polyline.
 *
 * The result is the sum of the lengths of all consecutive segments in
 * stored point order.
 *
 * Empty and singleton polylines have length zero.
 *
 * Uses the same `MetricScalar` policy as `segmentLength()`.
 *
 * Segment lengths are accumulated in stored order in `MetricScalar!T`
 * using Kahan-style compensated summation to reduce floating-point
 * accumulation error.
 *
 * Compensation improves mixed-scale sums but does not make the result
 * exact, correctly rounded, or independent of segment order. Each segment
 * length remains an ordinary floating-point metric computation.
 *
 * Non-finite segment lengths and accumulated overflow propagate according
 * to normal floating-point arithmetic. The accumulated result may therefore
 * be NaN or infinity.
 *
 * No allocation or point copying is performed.
 *
 * Complexity:
 *     O(n) time and O(1) auxiliary space for n stored points.
 */
MetricScalar!T polylineLength(T)(Polyline3View!T polyline)
    pure nothrow @safe @nogc
if (isGeoScalar!T)
{
    alias M = MetricScalar!T;

    M result = M(0);
    M correction = M(0);

    foreach (i; 0 .. polyline.segmentCount)
    {
        const M value =
            segmentLength(
                polyline.segment(i)
            );

        const M adjusted =
            value - correction;

        const M next =
            result + adjusted;

        /*
         * A non-finite next value covers:
         *
         * - a non-finite segment length;
         * - accumulated finite overflow;
         * - an already non-finite running result.
         *
         * Preserve that ordinary floating-point result and discard the
         * compensation state before continuing.
         */
        if (!isFinite(next))
        {
            result = next;
            correction = M(0);
            continue;
        }

        correction =
            (next - result) - adjusted;

        result = next;
    }

    return result;
}


/// Example summing consecutive three-dimensional segment lengths.
@safe unittest
{
    import geo3;

    alias P = Point3!double;

    P[3] points = [
        P(0.0, 0.0, 0.0),
        P(2.0, 3.0, 6.0),
        P(4.0, 6.0, 12.0)
    ];

    const polyline =
        Polyline3View!double(points[]);

    assert(
        polylineLength(polyline) ==
        14.0
    );
}


@safe unittest
{
    /*
     * Metric result policy remains shared with geo-d.
     */
    static assert(
        is(
            typeof(distance(
                Point3!int.init,
                Point3!int.init
            )) ==
            double
        )
    );

    static assert(
        is(
            typeof(distance(
                Point3!long.init,
                Point3!long.init
            )) ==
            double
        )
    );

    static assert(
        is(
            typeof(distance(
                Point3!float.init,
                Point3!float.init
            )) ==
            double
        )
    );

    static assert(
        is(
            typeof(distance(
                Point3!double.init,
                Point3!double.init
            )) ==
            double
        )
    );

    static assert(
        is(
            typeof(distance(
                Point3!real.init,
                Point3!real.init
            )) ==
            real
        )
    );


    /*
     * Signed metric differences preserve direction without signed
     * integer subtraction overflow.
     */
    assert(
        signedMetricDifference!int(5, 2) ==
        3.0
    );

    assert(
        signedMetricDifference!int(2, 5) ==
        -3.0
    );

    assert(
        signedMetricDifference!int(5, 5) ==
        0.0
    );

    assert(
        signedMetricDifference!long(
            long.max,
            long.max - 1
        ) ==
        1.0
    );

    assert(
        signedMetricDifference!long(
            long.max - 1,
            long.max
        ) ==
        -1.0
    );

    assert(
        signedMetricDifference!long(
            long.min + 1,
            long.min
        ) ==
        1.0
    );

    assert(
        signedMetricDifference!long(
            long.min,
            long.min + 1
        ) ==
        -1.0
    );

    /*
     * Crossing zero uses unsigned magnitudes rather than signed
     * subtraction.
     */
    assert(
        signedMetricDifference!long(1, -1) ==
        2.0
    );

    assert(
        signedMetricDifference!long(-1, 1) ==
        -2.0
    );

    assert(
        signedMetricDifference!long(
            long.max,
            long.min
        ) >
        0.0
    );

    assert(
        signedMetricDifference!long(
            long.min,
            long.max
        ) <
        0.0
    );


    /*
     * Basic three-dimensional metric behaviour.
     */
    alias P = Point3!double;

    const a =
        P(0.0, 0.0, 0.0);

    const b =
        P(2.0, 3.0, 6.0);

    assert(distance(a, b) == 7.0);
    assert(distance(b, a) == 7.0);
    assert(distance(a, a) == 0.0);


    /*
     * Small differences between very large long coordinates must survive.
     *
     * Each axis is tested independently so every 3D component is covered.
     */
    alias PL = Point3!long;

    assert(
        distance(
            PL(long.max, 0, 0),
            PL(long.max - 1, 0, 0)
        ) ==
        1.0
    );

    assert(
        distance(
            PL(0, long.max, 0),
            PL(0, long.max - 1, 0)
        ) ==
        1.0
    );

    assert(
        distance(
            PL(0, 0, long.max),
            PL(0, 0, long.max - 1)
        ) ==
        1.0
    );

    assert(
        distance(
            PL(long.min, 0, 0),
            PL(long.min + 1, 0, 0)
        ) ==
        1.0
    );

    assert(
        distance(
            PL(0, long.min, 0),
            PL(0, long.min + 1, 0)
        ) ==
        1.0
    );

    assert(
        distance(
            PL(0, 0, long.min),
            PL(0, 0, long.min + 1)
        ) ==
        1.0
    );


    /*
     * The complete signed-long span must not suffer signed subtraction
     * overflow. Precision after conversion to double is not claimed exact.
     */
    const extremeDistance =
        distance(
            PL(long.min, 0, 0),
            PL(long.max, 0, 0)
        );

    assert(extremeDistance > 0.0);
    assert(extremeDistance != double.infinity);


    /*
     * Signed int extremes receive the same overflow protection.
     */
    alias PI = Point3!int;

    assert(
        distance(
            PI(int.min, int.min, int.min),
            PI(int.max, int.max, int.max)
        ) >
        0.0
    );


    /*
     * Squared distance uses the same robust integral differencing policy.
     *
     * Cover each 3D axis independently so no component can accidentally
     * regress to convert-before-subtract behaviour.
     */
    assert(
        squaredDistance(
            PL(long.max, 0, 0),
            PL(long.max - 1, 0, 0)
        ) ==
        1.0
    );

    assert(
        squaredDistance(
            PL(0, long.max, 0),
            PL(0, long.max - 1, 0)
        ) ==
        1.0
    );

    assert(
        squaredDistance(
            PL(0, 0, long.max),
            PL(0, 0, long.max - 1)
        ) ==
        1.0
    );

    assert(
        squaredDistance(
            P(0.0, 0.0, 0.0),
            P(2.0, 3.0, 6.0)
        ) ==
        49.0
    );

    /*
     * squaredDistance deliberately permits ordinary floating-point square
     * overflow after the robust component differences have been obtained.
     */
    const hugeSquared =
        squaredDistance(
            P(-double.max, 0.0, 0.0),
            P(double.max, 0.0, 0.0)
        );

    assert(
        hugeSquared ==
        double.infinity
    );


    /*
     * Segment length is exactly the point-distance metric applied to the
     * stored endpoints, including 3D Z differences and degenerate segments.
     */
    alias S = Segment3!double;

    assert(
        segmentLength(
            S(
                P(0.0, 0.0, 0.0),
                P(2.0, 3.0, 6.0)
            )
        ) ==
        7.0
    );

    assert(
        segmentLength(
            S(
                P(4.0, -2.0, 9.0),
                P(4.0, -2.0, 9.0)
            )
        ) ==
        0.0
    );

    const longSegment =
        Segment3!long(
            PL(0, 0, long.max),
            PL(0, 0, long.max - 1)
        );

    assert(
        segmentLength(longSegment) ==
        1.0
    );


    /*
     * Polyline length is the sum of consecutive 3D segment lengths.
     */
    {
        alias PP = Point3!double;

        PP[3] points = [
            PP(0.0, 0.0, 0.0),
            PP(2.0, 3.0, 6.0),
            PP(4.0, 6.0, 12.0)
        ];

        auto polyline =
            Polyline3View!double(points[]);

        assert(polyline.segmentCount == 2);

        assert(
            polylineLength(polyline) ==
            14.0
        );
    }


    /*
     * Compensated accumulation preserves small segment lengths that ordinary
     * sequential addition can lose after the running total becomes large.
     *
     * Each block contributes exactly:
     *
     *     2 * 2^52 + 2
     *
     * and the complete expected result is itself exactly representable as
     * binary64.
     */
    {
        alias PP = Point3!double;

        enum size_t blocks = 256;
        enum double large = 0x1p52;
        enum double expected = 0x1p61 + 512.0;

        PP[1 + blocks * 4] points;

        size_t index;

        points[index++] =
            PP(0.0, 0.0, 0.0);

        foreach (_; 0 .. blocks)
        {
            points[index++] =
                PP(large, 0.0, 0.0);

            points[index++] =
                PP(0.0, 0.0, 0.0);

            points[index++] =
                PP(1.0, 0.0, 0.0);

            points[index++] =
                PP(0.0, 0.0, 0.0);
        }

        assert(index == points.length);

        const length =
            polylineLength(
                Polyline3View!double(points[])
            );

        assert(length == expected);
    }


    /*
     * Compensated accumulation retains normal floating-point non-finite
     * semantics.
     */
    {
        alias PP = Point3!double;

        /*
         * Exercise the third coordinate explicitly.
         */
        PP[2] infinitePoints = [
            PP(0.0, 0.0, 0.0),
            PP(0.0, 0.0, double.infinity)
        ];

        assert(
            polylineLength(
                Polyline3View!double(
                    infinitePoints[]
                )
            ) ==
            double.infinity
        );


        PP[2] nanPoints = [
            PP(0.0, 0.0, 0.0),
            PP(0.0, 0.0, double.nan)
        ];

        const nanLength =
            polylineLength(
                Polyline3View!double(
                    nanPoints[]
                )
            );

        assert(nanLength != nanLength);


        /*
         * Every individual segment length is finite, but two double.max
         * segments overflow the accumulated result. A later finite segment
         * must leave that infinity intact.
         */
        PP[4] overflowPoints = [
            PP(0.0, 0.0, 0.0),
            PP(double.max, 0.0, 0.0),
            PP(0.0, 0.0, 0.0),
            PP(1.0, 0.0, 0.0)
        ];

        assert(
            polylineLength(
                Polyline3View!double(
                    overflowPoints[]
                )
            ) ==
            double.infinity
        );
    }


    /*
     * Empty and singleton polylines have zero length.
     */
    {
        Point3!int[] emptyPoints;

        auto empty =
            Polyline3View!int(
                emptyPoints
            );

        assert(
            polylineLength(empty) ==
            0.0
        );

        Point3!int[1] singletonPoints = [
            Point3!int(7, -3, 11)
        ];

        auto singleton =
            Polyline3View!int(
                singletonPoints[]
            );

        assert(
            polylineLength(singleton) ==
            0.0
        );
    }


    /*
     * Polyline metric result types follow the shared MetricScalar policy.
     */
    static assert(
        is(
            typeof(
                polylineLength(
                    Polyline3View!int.init
                )
            ) ==
            double
        )
    );

    static assert(
        is(
            typeof(
                polylineLength(
                    Polyline3View!long.init
                )
            ) ==
            double
        )
    );

    static assert(
        is(
            typeof(
                polylineLength(
                    Polyline3View!float.init
                )
            ) ==
            double
        )
    );

    static assert(
        is(
            typeof(
                polylineLength(
                    Polyline3View!double.init
                )
            ) ==
            double
        )
    );

    static assert(
        is(
            typeof(
                polylineLength(
                    Polyline3View!real.init
                )
            ) ==
            real
        )
    );


    /*
     * Floating-point non-finite values retain normal metric semantics.
     */
    const infinitePoint =
        P(double.infinity, 0.0, 0.0);

    assert(
        distance(a, infinitePoint) ==
        double.infinity
    );

    const nanPoint =
        P(0.0, 0.0, double.nan);

    const nanDistance =
        distance(a, nanPoint);

    assert(nanDistance != nanDistance);
}
