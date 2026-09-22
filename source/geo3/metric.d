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

import std.math.exponential :
    ilogb,
    scalbn;

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


/*
 * Conservative domain for direct binary floating-point products used by the
 * ordinary-range point/segment metric fast path.
 *
 * Every non-zero component must lie in [2^-250, 2^250].
 *
 * With three dimensions, direct dot products remain below roughly 2^502 and
 * direct cross-product components below roughly 2^501, comfortably inside
 * the normal binary64 range.
 *
 * Values outside this domain use power-of-two normalization.
 */
private bool isDirectMetricComponent(M)(M value)
    pure nothrow @safe @nogc
{
    const M magnitude =
        value < M(0)
            ? -value
            : value;

    return
        magnitude == M(0) ||
        (
            magnitude >= M(0x1p-250) &&
            magnitude <= M(0x1p250)
        );
}


/*
 * True when every component needed by the direct 3D projection and
 * perpendicular-distance formulas lies in the conservative product range.
 */
private bool hasDirectMetricProductRange(M)(
    M dx,
    M dy,
    M dz,
    M rx,
    M ry,
    M rz
)
    pure nothrow @safe @nogc
{
    return
        isDirectMetricComponent(dx) &&
        isDirectMetricComponent(dy) &&
        isDirectMetricComponent(dz) &&
        isDirectMetricComponent(rx) &&
        isDirectMetricComponent(ry) &&
        isDirectMetricComponent(rz);
}


/*
 * Direct ordinary-range projection.
 *
 * Computes:
 *
 *     t = dot(r, d) / dot(d, d)
 *
 * Preconditions:
 *
 * - all components are finite;
 * - d is non-zero;
 * - hasDirectMetricProductRange(...) is true.
 */
private M directProjectionParameter(M)(
    M dx,
    M dy,
    M dz,
    M rx,
    M ry,
    M rz
)
    pure nothrow @safe @nogc
{
    const M numerator =
        rx * dx +
        ry * dy +
        rz * dz;

    const M denominator =
        dx * dx +
        dy * dy +
        dz * dz;

    return numerator / denominator;
}


/*
 * Direct ordinary-range perpendicular distance from r to the infinite line
 * through the origin with direction d.
 *
 * Computes:
 *
 *     ||d x r|| / ||d||
 *
 * Preconditions are identical to directProjectionParameter().
 */
private M directPerpendicularDistance(M)(
    M dx,
    M dy,
    M dz,
    M rx,
    M ry,
    M rz
)
    pure nothrow @safe @nogc
{
    const M cx =
        dy * rz -
        dz * ry;

    const M cy =
        dz * rx -
        dx * rz;

    const M cz =
        dx * ry -
        dy * rx;

    const M crossNorm =
        hypot(
            hypot(
                cx,
                cy
            ),
            cz
        );

    if (crossNorm == M(0))
        return M(0);

    const M directionNorm =
        hypot(
            hypot(
                dx,
                dy
            ),
            dz
        );

    return crossNorm / directionNorm;
}


/*
 * Projection parameter of r onto d.
 *
 * Computes:
 *
 *     t = dot(r, d) / dot(d, d)
 *
 * without directly forming potentially overflowing or underflowing dot
 * products outside the conservative direct-product domain.
 *
 * d and r are independently normalized by powers of two:
 *
 *     d = nd * 2^expD
 *     r = nr * 2^expR
 *
 * therefore:
 *
 *     t =
 *         dot(nr, nd) / dot(nd, nd)
 *         * 2^(expR - expD)
 *
 * Preconditions:
 *
 * - all components are finite;
 * - d is non-zero.
 */
private MetricScalar!T projectionParameter(T)(
    MetricScalar!T dx,
    MetricScalar!T dy,
    MetricScalar!T dz,
    MetricScalar!T rx,
    MetricScalar!T ry,
    MetricScalar!T rz
)
    pure nothrow @safe @nogc
if (isGeoScalar!T)
{
    alias M = MetricScalar!T;

    if (
        hasDirectMetricProductRange(
            dx,
            dy,
            dz,
            rx,
            ry,
            rz
        )
    )
    {
        return directProjectionParameter(
            dx,
            dy,
            dz,
            rx,
            ry,
            rz
        );
    }

    const M absDx = dx < M(0) ? -dx : dx;
    const M absDy = dy < M(0) ? -dy : dy;
    const M absDz = dz < M(0) ? -dz : dz;

    const M absRx = rx < M(0) ? -rx : rx;
    const M absRy = ry < M(0) ? -ry : ry;
    const M absRz = rz < M(0) ? -rz : rz;

    const M maxD =
        absDx > absDy
            ? (absDx > absDz ? absDx : absDz)
            : (absDy > absDz ? absDy : absDz);

    const M maxR =
        absRx > absRy
            ? (absRx > absRz ? absRx : absRz)
            : (absRy > absRz ? absRy : absRz);

    assert(maxD > M(0));

    if (maxR == M(0))
        return M(0);

    const int expD =
        ilogb(maxD);

    const int expR =
        ilogb(maxR);

    const M ndx =
        scalbn(dx, -expD);

    const M ndy =
        scalbn(dy, -expD);

    const M ndz =
        scalbn(dz, -expD);

    const M nrx =
        scalbn(rx, -expR);

    const M nry =
        scalbn(ry, -expR);

    const M nrz =
        scalbn(rz, -expR);

    const M numerator =
        nrx * ndx +
        nry * ndy +
        nrz * ndz;

    const M denominator =
        ndx * ndx +
        ndy * ndy +
        ndz * ndz;

    const M normalizedRatio =
        numerator /
        denominator;

    return scalbn(
        normalizedRatio,
        expR - expD
    );
}


/*
 * Perpendicular distance from r to the infinite line through the origin with
 * direction d.
 *
 * Computes:
 *
 *     ||d x r|| / ||d||
 *
 * without directly forming potentially overflowing products outside the
 * conservative direct-product domain.
 *
 * With independently normalized vectors
 *
 *     d = nd * 2^expD
 *     r = nr * 2^expR
 *
 * the direction scale cancels:
 *
 *     ||d x r|| / ||d||
 *
 *       = ||nd x nr|| / ||nd|| * 2^expR
 *
 * Preconditions:
 *
 * - all components are finite;
 * - d is non-zero.
 */
private MetricScalar!T perpendicularDistance(T)(
    MetricScalar!T dx,
    MetricScalar!T dy,
    MetricScalar!T dz,
    MetricScalar!T rx,
    MetricScalar!T ry,
    MetricScalar!T rz
)
    pure nothrow @safe @nogc
if (isGeoScalar!T)
{
    alias M = MetricScalar!T;

    if (
        hasDirectMetricProductRange(
            dx,
            dy,
            dz,
            rx,
            ry,
            rz
        )
    )
    {
        return directPerpendicularDistance(
            dx,
            dy,
            dz,
            rx,
            ry,
            rz
        );
    }

    const M absDx = dx < M(0) ? -dx : dx;
    const M absDy = dy < M(0) ? -dy : dy;
    const M absDz = dz < M(0) ? -dz : dz;

    const M absRx = rx < M(0) ? -rx : rx;
    const M absRy = ry < M(0) ? -ry : ry;
    const M absRz = rz < M(0) ? -rz : rz;

    const M maxD =
        absDx > absDy
            ? (absDx > absDz ? absDx : absDz)
            : (absDy > absDz ? absDy : absDz);

    const M maxR =
        absRx > absRy
            ? (absRx > absRz ? absRx : absRz)
            : (absRy > absRz ? absRy : absRz);

    assert(maxD > M(0));

    if (maxR == M(0))
        return M(0);

    const int expD =
        ilogb(maxD);

    const int expR =
        ilogb(maxR);

    const M ndx =
        scalbn(dx, -expD);

    const M ndy =
        scalbn(dy, -expD);

    const M ndz =
        scalbn(dz, -expD);

    const M nrx =
        scalbn(rx, -expR);

    const M nry =
        scalbn(ry, -expR);

    const M nrz =
        scalbn(rz, -expR);

    const M cx =
        ndy * nrz -
        ndz * nry;

    const M cy =
        ndz * nrx -
        ndx * nrz;

    const M cz =
        ndx * nry -
        ndy * nrx;

    const M crossNorm =
        hypot(
            hypot(
                cx,
                cy
            ),
            cz
        );

    /*
     * Preserve exact collinearity explicitly before scalbn().
     */
    if (crossNorm == M(0))
        return M(0);

    const M directionNorm =
        hypot(
            hypot(
                ndx,
                ndy
            ),
            ndz
        );

    const M normalizedDistance =
        crossNorm /
        directionNorm;

    return scalbn(
        normalizedDistance,
        expR
    );
}


// Internal regression coverage for the 3D point/segment metric helpers.
@safe unittest
{
    /*
     * Ordinary-range projection.
     *
     * d = (2, 0, 0)
     * r = (1, 3, 4)
     *
     * projection parameter = 1 / 2.
     */
    assert(
        directProjectionParameter(
            2.0,
            0.0,
            0.0,
            1.0,
            3.0,
            4.0
        ) ==
        0.5
    );

    assert(
        directPerpendicularDistance(
            2.0,
            0.0,
            0.0,
            1.0,
            3.0,
            4.0
        ) ==
        5.0
    );


    /*
     * Genuine 3D cross product: all three cross components are non-zero.
     *
     * d = (1, 2, 3)
     * r = (4, 5, 6)
     *
     * d x r = (-3, 6, -3)
     */
    const crossDistance =
        directPerpendicularDistance(
            1.0,
            2.0,
            3.0,
            4.0,
            5.0,
            6.0
        );

    assert(crossDistance > 1.9639);
    assert(crossDistance < 1.9640);


    /*
     * Large values force the scaled implementation.
     *
     * Naive dot products overflow, while the normalized projection remains:
     *
     *     t = 0.6
     */
    assert(
        !hasDirectMetricProductRange(
            1.0e300,
            1.0e300,
            1.0e300,
            5.0e299,
            6.0e299,
            7.0e299
        )
    );

    const hugeT =
        projectionParameter!double(
            1.0e300,
            1.0e300,
            1.0e300,
            5.0e299,
            6.0e299,
            7.0e299
        );

    assert(hugeT > 0.5999);
    assert(hugeT < 0.6001);


    /*
     * For the same vectors:
     *
     *     ||d x r|| / ||d|| = sqrt(2) * 1e299
     */
    const hugeDistance =
        perpendicularDistance!double(
            1.0e300,
            1.0e300,
            1.0e300,
            5.0e299,
            6.0e299,
            7.0e299
        );

    assert(hugeDistance > 1.413e299);
    assert(hugeDistance < 1.415e299);


    /*
     * Power-of-two-scaled exact collinearity must remain exact zero.
     *
     * This also exercises the explicit zero path before scalbn().
     */
    const collinearDistance =
        perpendicularDistance!double(
             0x1p800,
            -0x1.8p800,
             0x1p799,
             0x1p799,
            -0x1.8p799,
             0x1p798
        );

    assert(collinearDistance == 0.0);

    const collinearT =
        projectionParameter!double(
             0x1p800,
            -0x1.8p800,
             0x1p799,
             0x1p799,
            -0x1.8p799,
             0x1p798
        );

    assert(collinearT == 0.5);


    /*
     * A huge scale ratio may overflow the unclamped projection parameter.
     * This is intentional: the later public operation can still infer that
     * the nearest point is the endpoint at t >= 1.
     */
    const endpointT =
        projectionParameter!double(
            1.0e-300,
            0.0,
            0.0,
            1.0e300,
            0.0,
            0.0
        );

    assert(
        endpointT ==
        double.infinity
    );
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
 * Computes the Euclidean distance from a point to a segment.
 *
 * Returns false when:
 *
 * - an input coordinate is NaN or infinite; or
 * - a required metric difference cannot be represented finitely in
 *   `MetricScalar!T`.
 *
 * On failure, `result` is zero. A successful call produces a finite result.
 *
 * A degenerate segment is treated as its single endpoint.
 *
 * Integral coordinate differences are obtained before conversion to the
 * metric computation type, avoiding signed overflow and preserving small
 * differences between large integer coordinates.
 *
 * For an interior projection, the perpendicular distance is computed from
 * scaled three-dimensional direction and offset vectors. No rounded
 * nearest-point coordinate is constructed.
 *
 * This is a metric computation, not an exact topological predicate.
 *
 * No allocation is performed.
 *
 * Complexity:
 *     O(1) time and O(1) auxiliary space.
 */
bool tryPointSegmentDistance(T, R)(
    Segment3!T segment,
    Point3!T point,
    out R result
)
    pure nothrow @safe @nogc
if (
    isGeoScalar!T &&
    is(R == MetricScalar!T)
)
{
    alias M = MetricScalar!T;

    result = M(0);

    if (!point.isFinite || !segment.isFinite)
        return false;

    const M dx =
        signedMetricDifference(
            segment.b.x,
            segment.a.x
        );

    const M dy =
        signedMetricDifference(
            segment.b.y,
            segment.a.y
        );

    const M dz =
        signedMetricDifference(
            segment.b.z,
            segment.a.z
        );

    const M rx =
        signedMetricDifference(
            point.x,
            segment.a.x
        );

    const M ry =
        signedMetricDifference(
            point.y,
            segment.a.y
        );

    const M rz =
        signedMetricDifference(
            point.z,
            segment.a.z
        );

    if (
        !isFinite(dx) ||
        !isFinite(dy) ||
        !isFinite(dz) ||
        !isFinite(rx) ||
        !isFinite(ry) ||
        !isFinite(rz)
    )
    {
        return false;
    }

    /*
     * Degenerate segment.
     */
    if (
        dx == M(0) &&
        dy == M(0) &&
        dz == M(0)
    )
    {
        result =
            hypot(
                hypot(
                    rx,
                    ry
                ),
                rz
            );

        if (!isFinite(result))
        {
            result = M(0);
            return false;
        }

        return true;
    }

    const bool directMetricRange =
        hasDirectMetricProductRange(
            dx,
            dy,
            dz,
            rx,
            ry,
            rz
        );

    const M t =
        directMetricRange
            ? directProjectionParameter(
                dx,
                dy,
                dz,
                rx,
                ry,
                rz
            )
            : projectionParameter!T(
                dx,
                dy,
                dz,
                rx,
                ry,
                rz
            );

    /*
     * Infinite projection parameters still determine an endpoint.
     * NaN reaches the explicit finite check below.
     */
    if (t <= M(0))
    {
        result =
            hypot(
                hypot(
                    rx,
                    ry
                ),
                rz
            );

        if (!isFinite(result))
        {
            result = M(0);
            return false;
        }

        return true;
    }

    if (t >= M(1))
    {
        const M bx =
            signedMetricDifference(
                point.x,
                segment.b.x
            );

        const M by =
            signedMetricDifference(
                point.y,
                segment.b.y
            );

        const M bz =
            signedMetricDifference(
                point.z,
                segment.b.z
            );

        if (
            !isFinite(bx) ||
            !isFinite(by) ||
            !isFinite(bz)
        )
        {
            return false;
        }

        result =
            hypot(
                hypot(
                    bx,
                    by
                ),
                bz
            );

        if (!isFinite(result))
        {
            result = M(0);
            return false;
        }

        return true;
    }

    if (!isFinite(t))
        return false;

    result =
        directMetricRange
            ? directPerpendicularDistance(
                dx,
                dy,
                dz,
                rx,
                ry,
                rz
            )
            : perpendicularDistance!T(
                dx,
                dy,
                dz,
                rx,
                ry,
                rz
            );

    if (!isFinite(result))
    {
        result = M(0);
        return false;
    }

    return true;
}


/// Example computing a genuinely three-dimensional point-to-segment distance.
@safe unittest
{
    import geo3;

    alias P = Point3!double;
    alias S = Segment3!double;

    const segment =
        S(
            P(0.0, 0.0, 0.0),
            P(10.0, 0.0, 0.0)
        );

    double result;

    assert(
        tryPointSegmentDistance(
            segment,
            P(5.0, 3.0, 4.0),
            result
        )
    );

    assert(result == 5.0);
}


/**
 * Finds the nearest point on a segment to a point.
 *
 * The result uses `MetricScalar!T` because the nearest point of an integral
 * segment is not generally representable with integral coordinates.
 *
 * Returns false when:
 *
 * - an input coordinate is NaN or infinite; or
 * - a required metric difference cannot be represented finitely in the
 *   metric computation type.
 *
 * On failure, `result` remains `Point3!(MetricScalar!T).init`.
 * A successful result contains only finite coordinates.
 *
 * A degenerate segment returns its single endpoint converted to
 * `MetricScalar!T`. Endpoint projections are represented in the same way.
 * Consequently, long coordinates may be rounded when represented as double.
 *
 * Interior nearest points are constructed using floating-point metric
 * arithmetic. The constructed coordinates are not an exact topological
 * representation and must not be treated as an exact predicate result.
 *
 * No allocation is performed.
 *
 * Complexity:
 *     O(1) time and O(1) auxiliary space.
 */
bool tryNearestPoint(T, R)(
    Segment3!T segment,
    Point3!T point,
    out Point3!R result
)
    pure nothrow @safe @nogc
if (
    isGeoScalar!T &&
    is(R == MetricScalar!T)
)
{
    alias M = MetricScalar!T;

    if (!segment.isFinite || !point.isFinite)
        return false;

    const M ax = cast(M) segment.a.x;
    const M ay = cast(M) segment.a.y;
    const M az = cast(M) segment.a.z;

    const M bx = cast(M) segment.b.x;
    const M by = cast(M) segment.b.y;
    const M bz = cast(M) segment.b.z;

    const M dx =
        signedMetricDifference(
            segment.b.x,
            segment.a.x
        );

    const M dy =
        signedMetricDifference(
            segment.b.y,
            segment.a.y
        );

    const M dz =
        signedMetricDifference(
            segment.b.z,
            segment.a.z
        );

    if (
        !isFinite(dx) ||
        !isFinite(dy) ||
        !isFinite(dz)
    )
    {
        return false;
    }

    /*
     * Degenerate segment.
     */
    if (
        dx == M(0) &&
        dy == M(0) &&
        dz == M(0)
    )
    {
        result =
            Point3!M(
                ax,
                ay,
                az
            );

        return true;
    }

    const M rx =
        signedMetricDifference(
            point.x,
            segment.a.x
        );

    const M ry =
        signedMetricDifference(
            point.y,
            segment.a.y
        );

    const M rz =
        signedMetricDifference(
            point.z,
            segment.a.z
        );

    if (
        !isFinite(rx) ||
        !isFinite(ry) ||
        !isFinite(rz)
    )
    {
        return false;
    }

    const bool directMetricRange =
        hasDirectMetricProductRange(
            dx,
            dy,
            dz,
            rx,
            ry,
            rz
        );

    const M t =
        directMetricRange
            ? directProjectionParameter(
                dx,
                dy,
                dz,
                rx,
                ry,
                rz
            )
            : projectionParameter!T(
                dx,
                dy,
                dz,
                rx,
                ry,
                rz
            );

    /*
     * Positive overflow of t means that the projection lies beyond b.
     * Negative overflow means that it lies before a.
     */
    if (t <= M(0))
    {
        result =
            Point3!M(
                ax,
                ay,
                az
            );

        return true;
    }

    if (t >= M(1))
    {
        result =
            Point3!M(
                bx,
                by,
                bz
            );

        return true;
    }

    if (!isFinite(t))
        return false;

    /*
     * Interpolate from the nearer endpoint to reduce avoidable
     * cancellation when t is very close to zero or one.
     */
    if (t <= M(0.5))
    {
        result =
            Point3!M(
                ax + t * dx,
                ay + t * dy,
                az + t * dz
            );
    }
    else
    {
        const M fromB =
            t - M(1);

        result =
            Point3!M(
                bx + fromB * dx,
                by + fromB * dy,
                bz + fromB * dz
            );
    }

    if (!result.isFinite)
    {
        result = Point3!M.init;
        return false;
    }

    return true;
}


/// Example constructing the nearest point in the metric scalar type.
@safe unittest
{
    import geo3;

    const segment =
        Segment3!int(
            Point3!int(0, 0, 0),
            Point3!int(10, 0, 0)
        );

    Point3!double nearest;

    assert(
        tryNearestPoint(
            segment,
            Point3!int(3, 4, 5),
            nearest
        )
    );

    static assert(
        is(
            typeof(nearest) ==
            Point3!double
        )
    );

    assert(
        nearest ==
        Point3!double(
            3.0,
            0.0,
            0.0
        )
    );
}


// Public point/segment metric regression coverage.
@safe unittest
{
    alias P = Point3!double;
    alias S = Segment3!double;

    double d;
    Point3!double nearest;


    /*
     * Genuine 3D interior projection.
     */
    const diagonal =
        S(
            P(0.0, 0.0, 0.0),
            P(2.0, 2.0, 0.0)
        );

    assert(
        tryNearestPoint(
            diagonal,
            P(1.0, 1.0, 3.0),
            nearest
        )
    );

    assert(
        nearest ==
        P(1.0, 1.0, 0.0)
    );

    assert(
        tryPointSegmentDistance(
            diagonal,
            P(1.0, 1.0, 3.0),
            d
        )
    );

    assert(d > 2.999999999999999);
    assert(d < 3.000000000000001);


    /*
     * Projection before A.
     */
    const axis =
        S(
            P(0.0, 0.0, 0.0),
            P(10.0, 0.0, 0.0)
        );

    assert(
        tryNearestPoint(
            axis,
            P(-3.0, 4.0, 12.0),
            nearest
        )
    );

    assert(
        nearest ==
        P(0.0, 0.0, 0.0)
    );

    assert(
        tryPointSegmentDistance(
            axis,
            P(-3.0, 4.0, 12.0),
            d
        )
    );

    assert(d == 13.0);


    /*
     * Projection beyond B.
     */
    assert(
        tryNearestPoint(
            axis,
            P(13.0, 4.0, 12.0),
            nearest
        )
    );

    assert(
        nearest ==
        P(10.0, 0.0, 0.0)
    );

    assert(
        tryPointSegmentDistance(
            axis,
            P(13.0, 4.0, 12.0),
            d
        )
    );

    assert(d == 13.0);


    /*
     * Degenerate segment behaves as a point.
     */
    const degenerate =
        S(
            P(1.0, 2.0, 3.0),
            P(1.0, 2.0, 3.0)
        );

    assert(
        tryNearestPoint(
            degenerate,
            P(4.0, 6.0, 15.0),
            nearest
        )
    );

    assert(
        nearest ==
        P(1.0, 2.0, 3.0)
    );

    assert(
        tryPointSegmentDistance(
            degenerate,
            P(4.0, 6.0, 15.0),
            d
        )
    );

    assert(d == 13.0);


    /*
     * Large products use the scaled path.
     */
    const huge =
        S(
            P(0.0, 0.0, 0.0),
            P(1.0e300, 1.0e300, 1.0e300)
        );

    assert(
        tryNearestPoint(
            huge,
            P(
                5.0e299,
                6.0e299,
                7.0e299
            ),
            nearest
        )
    );

    assert(nearest.x > 5.99e299);
    assert(nearest.x < 6.01e299);
    assert(nearest.y > 5.99e299);
    assert(nearest.y < 6.01e299);
    assert(nearest.z > 5.99e299);
    assert(nearest.z < 6.01e299);

    assert(
        tryPointSegmentDistance(
            huge,
            P(
                5.0e299,
                6.0e299,
                7.0e299
            ),
            d
        )
    );

    assert(d > 1.413e299);
    assert(d < 1.415e299);


    /*
     * Very different scales may overflow the unclamped projection
     * parameter, while endpoint clamping remains well-defined.
     */
    const tiny =
        S(
            P(0.0, 0.0, 0.0),
            P(1.0e-300, 0.0, 0.0)
        );

    assert(
        tryNearestPoint(
            tiny,
            P(1.0e300, 0.0, 0.0),
            nearest
        )
    );

    assert(
        nearest ==
        P(1.0e-300, 0.0, 0.0)
    );


    /*
     * Non-finite geometry is rejected transactionally.
     */
    nearest =
        P(99.0, 99.0, 99.0);

    assert(
        !tryNearestPoint(
            axis,
            P(
                0.0,
                0.0,
                double.nan
            ),
            nearest
        )
    );

    assert(
        nearest ==
        Point3!double.init
    );

    d = 123.0;

    assert(
        !tryPointSegmentDistance(
            axis,
            P(
                0.0,
                double.infinity,
                0.0
            ),
            d
        )
    );

    assert(d == 0.0);


    /*
     * Finite inputs whose component difference exceeds MetricScalar range
     * are conservatively rejected.
     */
    const extreme =
        S(
            P(
                -double.max,
                0.0,
                0.0
            ),
            P(
                double.max,
                0.0,
                0.0
            )
        );

    assert(
        !tryNearestPoint(
            extreme,
            P(0.0, 1.0, 1.0),
            nearest
        )
    );

    assert(
        nearest ==
        Point3!double.init
    );

    d = 123.0;

    assert(
        !tryPointSegmentDistance(
            extreme,
            P(0.0, 1.0, 1.0),
            d
        )
    );

    assert(d == 0.0);


    /*
     * Large integral coordinates retain small differences before conversion.
     */
    {
        alias PL = Point3!long;
        alias SL = Segment3!long;

        double integralDistance;
        Point3!double integralNearest;

        const long base =
            long.max - 8192;

        assert(
            tryNearestPoint(
                SL(
                    PL(base, 0, 0),
                    PL(base + 8192, 0, 0)
                ),
                PL(
                    base + 4096,
                    100,
                    -200
                ),
                integralNearest
            )
        );

        assert(integralNearest.y == 0.0);
        assert(integralNearest.z == 0.0);

        assert(
            tryPointSegmentDistance(
                SL(
                    PL(long.min, 0, 0),
                    PL(long.max, 0, 0)
                ),
                PL(0, 3, 4),
                integralDistance
            )
        );

        assert(integralDistance == 5.0);
    }


    /*
     * Metric result policy.
     */
    static assert(
        is(
            typeof({
                double value;
                tryPointSegmentDistance(
                    Segment3!int.init,
                    Point3!int.init,
                    value
                );
                return value;
            }()) ==
            double
        )
    );

    static assert(
        is(
            typeof({
                real value;
                tryPointSegmentDistance(
                    Segment3!real.init,
                    Point3!real.init,
                    value
                );
                return value;
            }()) ==
            real
        )
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
