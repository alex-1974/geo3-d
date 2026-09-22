/**
 * Axis-aligned bounds of three-dimensional geometry representations.
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
module geo3.bounding_box;

import geo3.bounds :
    Bounds3;

import geo3.linear_ring_view :
    LinearRing3View;

import geo3.point :
    Point3;

import geo3.polyline_view :
    Polyline3View;

import geo3.scalar :
    isGeoScalar;

import geo3.segment :
    Segment3;

import std.traits :
    isFloatingPoint;


/*
 * Incremental accumulator used by the public tryBounds overloads.
 *
 * The accumulator is intentionally private. It permits complete geometry
 * traversal before constructing the public Bounds3 result, preserving the
 * transactional out-parameter semantics when a later point contains NaN.
 */
private struct BoundsAccumulator(T)
if (isGeoScalar!T)
{
    private bool hasValue;

    private T minX;
    private T minY;
    private T minZ;

    private T maxX;
    private T maxY;
    private T maxZ;


    bool tryAdd(Point3!T point)
        pure nothrow @safe @nogc
    {
        static if (isFloatingPoint!T)
        {
            if (
                point.x != point.x ||
                point.y != point.y ||
                point.z != point.z
            )
            {
                return false;
            }
        }

        if (!hasValue)
        {
            minX = point.x;
            minY = point.y;
            minZ = point.z;

            maxX = point.x;
            maxY = point.y;
            maxZ = point.z;

            hasValue = true;

            return true;
        }

        if (point.x < minX)
            minX = point.x;

        if (point.y < minY)
            minY = point.y;

        if (point.z < minZ)
            minZ = point.z;

        if (point.x > maxX)
            maxX = point.x;

        if (point.y > maxY)
            maxY = point.y;

        if (point.z > maxZ)
            maxZ = point.z;

        return true;
    }


    bool finish(out Bounds3!T result) const
        pure nothrow @safe @nogc
    {
        if (!hasValue)
            return true;

        return Bounds3!T.tryFromMinMax(
            Point3!T(
                minX,
                minY,
                minZ
            ),
            Point3!T(
                maxX,
                maxY,
                maxZ
            ),
            result
        );
    }
}


/**
 * Computes the closed axis-aligned bounds of a segment.
 *
 * Both stored endpoints participate in the result. Endpoint order does not
 * affect the geometric bounds.
 *
 * A degenerate segment produces a degenerate non-empty bounds.
 *
 * Returns false when either endpoint contains NaN. Infinity is permitted.
 *
 * On failure, `result` is `Bounds3!T.init`.
 *
 * No allocation is performed.
 *
 * Complexity:
 *
 *     time  O(1)
 *     space O(1)
 */
bool tryBounds(T)(
    Segment3!T segment,
    out Bounds3!T result
)
    pure nothrow @safe @nogc
if (isGeoScalar!T)
{
    BoundsAccumulator!T accumulated;

    if (!accumulated.tryAdd(segment.a))
        return false;

    if (!accumulated.tryAdd(segment.b))
        return false;

    return accumulated.finish(result);
}


/**
 * Computes the closed axis-aligned bounds of a polyline.
 *
 * Every stored point participates in the result.
 *
 * An empty polyline succeeds with `Bounds3!T.init`.
 * A singleton polyline produces a degenerate non-empty bounds.
 *
 * Returns false when any stored point contains NaN. Infinity is permitted.
 *
 * Failure is transactional: on failure, `result` is `Bounds3!T.init`.
 *
 * No allocation is performed.
 *
 * Complexity for `n` stored points:
 *
 *     time  O(n)
 *     space O(1)
 */
bool tryBounds(T)(
    scope Polyline3View!T polyline,
    out Bounds3!T result
)
    pure nothrow @safe @nogc
if (isGeoScalar!T)
{
    BoundsAccumulator!T accumulated;

    foreach (i; 0 .. polyline.length)
    {
        if (!accumulated.tryAdd(polyline[i]))
            return false;
    }

    return accumulated.finish(result);
}


/**
 * Computes the closed axis-aligned bounds of a linear ring.
 *
 * Every stored vertex participates in the result. The implicit closing
 * segment requires no separate treatment because both of its endpoints are
 * already stored vertices.
 *
 * Ring topology and planarity are not validated. Empty, degenerate,
 * non-planar, and otherwise geometrically invalid rings may still have
 * well-defined axis-aligned bounds.
 *
 * An empty ring succeeds with `Bounds3!T.init`.
 *
 * Returns false when any stored vertex contains NaN. Infinity is permitted.
 *
 * Failure is transactional: on failure, `result` is `Bounds3!T.init`.
 *
 * No allocation is performed.
 *
 * Complexity for `n` stored vertices:
 *
 *     time  O(n)
 *     space O(1)
 */
bool tryBounds(T)(
    scope LinearRing3View!T ring,
    out Bounds3!T result
)
    pure nothrow @safe @nogc
if (isGeoScalar!T)
{
    BoundsAccumulator!T accumulated;

    foreach (i; 0 .. ring.length)
    {
        if (!accumulated.tryAdd(ring[i]))
            return false;
    }

    return accumulated.finish(result);
}


/// Example computing bounds of a non-owning three-dimensional polyline view.
@safe unittest
{
    import geo3;

    alias P = Point3!int;

    P[4] points = [
        P(3,  4,  7),
        P(-2, 9, -1),
        P(8, -5, 12),
        P(1,  6,  3)
    ];

    const polyline =
        Polyline3View!int(points[]);

    Bounds3!int bounds;

    assert(
        tryBounds(
            polyline,
            bounds
        )
    );

    assert(
        bounds.min ==
        P(-2, -5, -1)
    );

    assert(
        bounds.max ==
        P(8, 9, 12)
    );
}


// Regression coverage for the complete scalar and representation family.
@safe unittest
{
    import std.meta :
        AliasSeq;

    static foreach (T; AliasSeq!(int, long, float, double, real))
    {
        {
            alias P = Point3!T;
        alias B = Bounds3!T;
        alias S = Segment3!T;
        alias L = Polyline3View!T;
        alias R = LinearRing3View!T;


        /*
         * Ordinary and reversed segments have identical bounds.
         */
        const S segment =
            S(
                P(T(4), T(-2), T(7)),
                P(T(-3), T(8), T(-5))
            );

        const S reversedSegment =
            S(
                segment.b,
                segment.a
            );

        B segmentBounds;
        B reversedSegmentBounds;

        assert(
            tryBounds(
                segment,
                segmentBounds
            )
        );

        assert(
            tryBounds(
                reversedSegment,
                reversedSegmentBounds
            )
        );

        assert(segmentBounds == reversedSegmentBounds);

        assert(
            segmentBounds.min ==
            P(T(-3), T(-2), T(-5))
        );

        assert(
            segmentBounds.max ==
            P(T(4), T(8), T(7))
        );


        /*
         * Degenerate segment.
         */
        B degenerateBounds;

        assert(
            tryBounds(
                S(
                    P(T(2), T(3), T(4)),
                    P(T(2), T(3), T(4))
                ),
                degenerateBounds
            )
        );

        assert(!degenerateBounds.empty);

        assert(
            degenerateBounds.min ==
            P(T(2), T(3), T(4))
        );

        assert(
            degenerateBounds.max ==
            P(T(2), T(3), T(4))
        );


        /*
         * Empty polyline has empty bounds.
         */
        P[] noPoints;

        auto emptyPolyline =
            L(noPoints);

        B emptyPolylineBounds;

        assert(
            tryBounds(
                emptyPolyline,
                emptyPolylineBounds
            )
        );

        assert(
            emptyPolylineBounds ==
            B.init
        );


        /*
         * Singleton polyline has degenerate bounds.
         */
        P[1] singletonPoints = [
            P(T(5), T(-7), T(11))
        ];

        auto singletonPolyline =
            L(singletonPoints[]);

        B singletonBounds;

        assert(
            tryBounds(
                singletonPolyline,
                singletonBounds
            )
        );

        assert(
            singletonBounds.min ==
            singletonPoints[0]
        );

        assert(
            singletonBounds.max ==
            singletonPoints[0]
        );


        /*
         * Ordinary polyline. Every axis reaches its extrema at potentially
         * different stored points.
         */
        P[4] polylinePoints = [
            P(T(3),  T(4),  T(7)),
            P(T(-2), T(9),  T(-6)),
            P(T(8),  T(-5), T(2)),
            P(T(1),  T(6),  T(12))
        ];

        auto polyline =
            L(polylinePoints[]);

        B polylineBounds;

        assert(
            tryBounds(
                polyline,
                polylineBounds
            )
        );

        assert(
            polylineBounds.min ==
            P(T(-2), T(-5), T(-6))
        );

        assert(
            polylineBounds.max ==
            P(T(8), T(9), T(12))
        );

        foreach (point; polylinePoints)
            assert(polylineBounds.contains(point));


        /*
         * Reversing stored point order does not change bounds.
         */
        P[4] reversedPolylinePoints = [
            polylinePoints[3],
            polylinePoints[2],
            polylinePoints[1],
            polylinePoints[0]
        ];

        B reversedPolylineBounds;

        assert(
            tryBounds(
                L(reversedPolylinePoints[]),
                reversedPolylineBounds
            )
        );

        assert(
            reversedPolylineBounds ==
            polylineBounds
        );


        /*
         * Empty ring has empty bounds.
         */
        B emptyRingBounds;

        assert(
            tryBounds(
                R(noPoints),
                emptyRingBounds
            )
        );

        assert(
            emptyRingBounds ==
            B.init
        );


        /*
         * Ring bounds depend only on stored vertices.
         *
         * The points deliberately form a non-planar cyclic sequence.
         * Planarity is not part of the representation or bounds contract.
         */
        P[4] ringPoints = [
            P(T(-4), T(1),  T(3)),
            P(T(6),  T(2),  T(-8)),
            P(T(3),  T(10), T(12)),
            P(T(-1), T(7),  T(5))
        ];

        auto ring =
            R(ringPoints[]);

        B ringBounds;

        assert(
            tryBounds(
                ring,
                ringBounds
            )
        );

        assert(
            ringBounds.min ==
            P(T(-4), T(1), T(-8))
        );

        assert(
            ringBounds.max ==
            P(T(6), T(10), T(12))
        );

        foreach (point; ringPoints)
            assert(ringBounds.contains(point));


        /*
         * Reversing cyclic traversal does not change bounds.
         */
        P[4] reversedRingPoints = [
            ringPoints[3],
            ringPoints[2],
            ringPoints[1],
            ringPoints[0]
        ];

        B reversedRingBounds;

        assert(
            tryBounds(
                R(reversedRingPoints[]),
                reversedRingBounds
            )
        );

        assert(
            reversedRingBounds ==
            ringBounds
        );


        /*
         * Integral extrema require no arithmetic and are preserved exactly.
         */
        static if (is(T == int) || is(T == long))
        {
            B extremeBounds;

            assert(
                tryBounds(
                    S(
                        P(T.min, T.max, T.min),
                        P(T.max, T.min, T.max)
                    ),
                    extremeBounds
                )
            );

            assert(
                extremeBounds.min ==
                P(T.min, T.min, T.min)
            );

            assert(
                extremeBounds.max ==
                P(T.max, T.max, T.max)
            );
        }


        static if (isFloatingPoint!T)
        {
            const T infinity =
                cast(T) double.infinity;

            /*
             * Infinity is valid and contributes normally, including on Z.
             */
            P[2] infinitePoints = [
                P(T(1), T(2), -infinity),
                P(T(3), T(4),  infinity)
            ];

            B infiniteBounds;

            assert(
                tryBounds(
                    L(infinitePoints[]),
                    infiniteBounds
                )
            );

            assert(!infiniteBounds.empty);
            assert(!infiniteBounds.isFinite);

            assert(
                infiniteBounds.min.z ==
                -infinity
            );

            assert(
                infiniteBounds.max.z ==
                infinity
            );


            /*
             * NaN after valid points must fail transactionally.
             *
             * The out parameter begins as a non-empty sentinel. Failure must
             * leave Bounds3.init rather than the sentinel or partial bounds.
             */
            P[3] nanPoints = [
                P(T(1), T(2), T(3)),
                P(T(4), T(5), T(6)),
                P(T(7), T(8), cast(T) double.nan)
            ];

            B failedBounds;

            assert(
                B.tryFromPoint(
                    P(T(99), T(99), T(99)),
                    failedBounds
                )
            );

            assert(!failedBounds.empty);

            assert(
                !tryBounds(
                    L(nanPoints[]),
                    failedBounds
                )
            );

            assert(
                failedBounds ==
                B.init
            );


            /*
             * Ring failure is equally transactional.
             */
            P[3] nanRingPoints = [
                P(T(20), T(20), T(20)),
                P(T(21), T(20), T(21)),
                P(T(20), T(21), cast(T) double.nan)
            ];

            assert(
                B.tryFromPoint(
                    P(T(99), T(99), T(99)),
                    failedBounds
                )
            );

            assert(
                !tryBounds(
                    R(nanRingPoints[]),
                    failedBounds
                )
            );

            assert(
                failedBounds ==
                B.init
            );


            /*
             * Segment NaN rejection.
             */
            assert(
                B.tryFromPoint(
                    P(T(99), T(99), T(99)),
                    failedBounds
                )
            );

            assert(
                !tryBounds(
                    S(
                        P(T(0), T(0), T(0)),
                        P(T(1), T(2), cast(T) double.nan)
                    ),
                    failedBounds
                )
            );

            assert(
                failedBounds ==
                B.init
            );
        }
        }
    }
}
