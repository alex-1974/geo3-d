/**
 * Certified first-stage binary64 filter for three-dimensional orientation.
 *
 * INTERNAL IMPLEMENTATION MODULE.
 *
 * The filter evaluates the classical orient3d determinant layout relative
 * to d and uses Shewchuk's first-stage error bound.
 *
 * A returned negative/coplanar/positive state is certified.
 *
 * `uncertain` is an internal control-flow state only. It means that a later
 * exact backend must determine the public result.
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
module geo3.internal.orientation_filter;

import geo3.internal.binary64_rounding :
    roundedAdd,
    roundedMul,
    roundedSub;

import geo3.point :
    Point3;

import std.bitmanip :
    DoubleRep;


/**
 * Internal result of the certified binary64 orient3d filter.
 *
 * negative/coplanar/positive use geo3-d's public determinant convention:
 *
 *     sign(det(b-a, c-a, d-a))
 */
enum OrientationFilterResult : byte
{
    negative  = -1,
    coplanar  =  0,
    positive  =  1,
    uncertain =  2,
}


/*
 * Shewchuk orient3d first-stage error coefficient:
 *
 *     epsilon = 2^-53
 *
 *     o3derrboundA
 *         = (7 + 56 * epsilon) * epsilon
 *
 * Exact binary64 representation:
 *
 *     0x1.c000000000007p-51
 */
private enum double o3dErrboundA =
    0x1.c000000000007p-51;


/*
 * Binary64 finiteness test using the IEEE-754 exponent field directly.
 *
 * Raw exponent 2047 denotes infinity or NaN.
 */
pragma(inline, true)
private bool finiteBits(double value)
    pure nothrow @safe @nogc
{
    DoubleRep representation;

    representation.value =
        value;

    return
        representation.exponent !=
        0x7ff;
}


/*
 * True exactly for binary64 normals and signed zeros.
 *
 * Subnormal values, infinity, and NaN return false.
 */
pragma(inline, true)
private bool normalOrZeroBits(double value)
    pure nothrow @safe @nogc
{
    DoubleRep representation;

    representation.value =
        value;


    const uint exponent =
        representation.exponent;


    if (exponent == 0)
    {
        return
            representation.fraction ==
            0;
    }


    return
        exponent !=
        0x7ff;
}


/*
 * Valid multiplication result for the certified first-stage filter.
 *
 * Normal products are accepted.
 *
 * Exact zero is accepted only when at least one operand is mathematically
 * zero. If both operands are non-zero but the rounded product is zero,
 * multiplication underflow occurred and the filter must fall back.
 *
 * Subnormal and non-finite products are rejected.
 */
pragma(inline, true)
private bool validProductBits(
    double lhs,
    double rhs,
    double product
)
    pure nothrow @safe @nogc
{
    DoubleRep representation;

    representation.value =
        product;


    const uint exponent =
        representation.exponent;


    if (exponent == 0)
    {
        if (
            representation.fraction !=
            0
        )
        {
            return false;
        }


        return
            lhs == 0.0 ||
            rhs == 0.0;
    }


    return
        exponent !=
        0x7ff;
}


/*
 * The certified error bound must be finite, normal, and non-zero.
 */
pragma(inline, true)
private bool validErrorBoundBits(double value)
    pure nothrow @safe @nogc
{
    DoubleRep representation;

    representation.value =
        value;


    const uint exponent =
        representation.exponent;


    return
        exponent != 0 &&
        exponent != 0x7ff;
}


/*
 * Absolute value without introducing an external math operation.
 *
 * Inputs reaching this helper are finite.
 */
pragma(inline, true)
private double absolute(double value)
    pure nothrow @safe @nogc
{
    return value < 0.0
        ? -value
        : value;
}


/*
 * Finite-coordinate guard for the defensive filter entry point.
 */
pragma(inline, true)
private bool finitePoint(Point3!double point)
    pure nothrow @safe @nogc
{
    return
        finiteBits(point.x) &&
        finiteBits(point.y) &&
        finiteBits(point.z);
}


/**
 * Defensive first-stage binary64 Orientation3 filter.
 *
 * This entry point preserves the historical internal behaviour:
 * non-finite coordinates are rejected conservatively as `uncertain`.
 *
 * Public `orientation(Point3!double, ...)` already has finite coordinates as
 * a precondition and therefore calls `orientationFilterFinite` directly to
 * avoid repeating these checks on the hot path.
 */
OrientationFilterResult orientationFilter(
    Point3!double a,
    Point3!double b,
    Point3!double c,
    Point3!double d
)
    pure nothrow @safe @nogc
{
    if (
        !finitePoint(a) ||
        !finitePoint(b) ||
        !finitePoint(c) ||
        !finitePoint(d)
    )
    {
        return
            OrientationFilterResult.uncertain;
    }


    return orientationFilterFinite(
        a,
        b,
        c,
        d
    );
}


/**
 * Certified first-stage binary64 Orientation3 filter for finite coordinates.
 *
 * Preconditions:
 *
 *     every input coordinate is finite.
 *
 * The public geo3-d convention is:
 *
 *     sign(det(b-a, c-a, d-a))
 *
 * Internally this function uses the classical orient3d layout relative to d:
 *
 *     det(a-d, b-d, c-d)
 *
 * These determinants have opposite sign, so every certified non-zero result
 * is sign-normalized before return.
 *
 * All intermediate classification uses direct IEEE-754 binary64 bit
 * inspection. Subnormal values, arithmetic overflow, multiplication
 * underflow, and values too close to the error bound conservatively return
 * `uncertain`.
 *
 * A `coplanar` result is returned only when the permanent is exactly zero
 * after all potentially underflowing multiplications have been rejected.
 */
OrientationFilterResult orientationFilterFinite(
    Point3!double a,
    Point3!double b,
    Point3!double c,
    Point3!double d
)
    pure nothrow @safe @nogc
{
    /*
     * Classical orient3d differences relative to d.
     */
    const double adx =
        roundedSub(
            a.x,
            d.x
        );

    const double bdx =
        roundedSub(
            b.x,
            d.x
        );

    const double cdx =
        roundedSub(
            c.x,
            d.x
        );

    const double ady =
        roundedSub(
            a.y,
            d.y
        );

    const double bdy =
        roundedSub(
            b.y,
            d.y
        );

    const double cdy =
        roundedSub(
            c.y,
            d.y
        );

    const double adz =
        roundedSub(
            a.z,
            d.z
        );

    const double bdz =
        roundedSub(
            b.z,
            d.z
        );

    const double cdz =
        roundedSub(
            c.z,
            d.z
        );


    if (
        !normalOrZeroBits(adx) ||
        !normalOrZeroBits(bdx) ||
        !normalOrZeroBits(cdx) ||
        !normalOrZeroBits(ady) ||
        !normalOrZeroBits(bdy) ||
        !normalOrZeroBits(cdy) ||
        !normalOrZeroBits(adz) ||
        !normalOrZeroBits(bdz) ||
        !normalOrZeroBits(cdz)
    )
    {
        return
            OrientationFilterResult.uncertain;
    }


    const double bdxcdy =
        roundedMul(
            bdx,
            cdy
        );

    const double cdxbdy =
        roundedMul(
            cdx,
            bdy
        );

    const double cdxady =
        roundedMul(
            cdx,
            ady
        );

    const double adxcdy =
        roundedMul(
            adx,
            cdy
        );

    const double adxbdy =
        roundedMul(
            adx,
            bdy
        );

    const double bdxady =
        roundedMul(
            bdx,
            ady
        );


    if (
        !validProductBits(
            bdx,
            cdy,
            bdxcdy
        ) ||
        !validProductBits(
            cdx,
            bdy,
            cdxbdy
        ) ||
        !validProductBits(
            cdx,
            ady,
            cdxady
        ) ||
        !validProductBits(
            adx,
            cdy,
            adxcdy
        ) ||
        !validProductBits(
            adx,
            bdy,
            adxbdy
        ) ||
        !validProductBits(
            bdx,
            ady,
            bdxady
        )
    )
    {
        return
            OrientationFilterResult.uncertain;
    }


    const double bc =
        roundedSub(
            bdxcdy,
            cdxbdy
        );

    const double ca =
        roundedSub(
            cdxady,
            adxcdy
        );

    const double ab =
        roundedSub(
            adxbdy,
            bdxady
        );


    if (
        !normalOrZeroBits(bc) ||
        !normalOrZeroBits(ca) ||
        !normalOrZeroBits(ab)
    )
    {
        return
            OrientationFilterResult.uncertain;
    }


    const double termA =
        roundedMul(
            adz,
            bc
        );

    const double termB =
        roundedMul(
            bdz,
            ca
        );

    const double termC =
        roundedMul(
            cdz,
            ab
        );


    if (
        !validProductBits(
            adz,
            bc,
            termA
        ) ||
        !validProductBits(
            bdz,
            ca,
            termB
        ) ||
        !validProductBits(
            cdz,
            ab,
            termC
        )
    )
    {
        return
            OrientationFilterResult.uncertain;
    }


    const double detAB =
        roundedAdd(
            termA,
            termB
        );

    if (!normalOrZeroBits(detAB))
    {
        return
            OrientationFilterResult.uncertain;
    }


    const double det =
        roundedAdd(
            detAB,
            termC
        );

    if (!normalOrZeroBits(det))
    {
        return
            OrientationFilterResult.uncertain;
    }


    const double permanentPairA =
        roundedAdd(
            absolute(bdxcdy),
            absolute(cdxbdy)
        );

    const double permanentPairB =
        roundedAdd(
            absolute(cdxady),
            absolute(adxcdy)
        );

    const double permanentPairC =
        roundedAdd(
            absolute(adxbdy),
            absolute(bdxady)
        );


    if (
        !normalOrZeroBits(permanentPairA) ||
        !normalOrZeroBits(permanentPairB) ||
        !normalOrZeroBits(permanentPairC)
    )
    {
        return
            OrientationFilterResult.uncertain;
    }


    const double absAdz =
        absolute(adz);

    const double absBdz =
        absolute(bdz);

    const double absCdz =
        absolute(cdz);


    const double permanentA =
        roundedMul(
            permanentPairA,
            absAdz
        );

    const double permanentB =
        roundedMul(
            permanentPairB,
            absBdz
        );

    const double permanentC =
        roundedMul(
            permanentPairC,
            absCdz
        );


    if (
        !validProductBits(
            permanentPairA,
            absAdz,
            permanentA
        ) ||
        !validProductBits(
            permanentPairB,
            absBdz,
            permanentB
        ) ||
        !validProductBits(
            permanentPairC,
            absCdz,
            permanentC
        )
    )
    {
        return
            OrientationFilterResult.uncertain;
    }


    const double permanentAB =
        roundedAdd(
            permanentA,
            permanentB
        );

    if (!normalOrZeroBits(permanentAB))
    {
        return
            OrientationFilterResult.uncertain;
    }


    const double permanent =
        roundedAdd(
            permanentAB,
            permanentC
        );

    if (!normalOrZeroBits(permanent))
    {
        return
            OrientationFilterResult.uncertain;
    }


    if (permanent == 0.0)
    {
        return
            OrientationFilterResult.coplanar;
    }


    const double errbound =
        roundedMul(
            o3dErrboundA,
            permanent
        );


    if (!validErrorBoundBits(errbound))
    {
        return
            OrientationFilterResult.uncertain;
    }


    /*
     * `det` uses the classical orient3d convention.
     *
     * geo3-d uses the opposite public convention, hence the inversion.
     */
    if (det > errbound)
    {
        return
            OrientationFilterResult.negative;
    }


    if (-det > errbound)
    {
        return
            OrientationFilterResult.positive;
    }


    return
        OrientationFilterResult.uncertain;
}



static assert(
    cast(byte) OrientationFilterResult.negative == -1
);

static assert(
    cast(byte) OrientationFilterResult.coplanar == 0
);

static assert(
    cast(byte) OrientationFilterResult.positive == 1
);

static assert(
    cast(byte) OrientationFilterResult.uncertain == 2
);


version(unittest)
{
    import std.bigint :
        BigInt;

    import std.bitmanip :
        DoubleRep;

    import std.math.traits :
        isFinite;


    private alias P =
        Point3!double;


    /*
     * Decode one finite binary64 value exactly as an integer multiple
     * of 2^-1074.
     */
    private BigInt oracleBinary64Integer(double value)
        @safe
    {
        assert(isFinite(value));

        DoubleRep representation;

        representation.value =
            value;

        const ulong fraction =
            representation.fraction;

        const uint rawExponent =
            representation.exponent;

        ulong mantissa;
        uint shift;

        if (rawExponent == 0)
        {
            mantissa =
                fraction;

            shift = 0;
        }
        else
        {
            mantissa =
                (1UL << 52) |
                fraction;

            shift =
                rawExponent - 1;
        }

        if (mantissa == 0)
            return BigInt(0);

        BigInt result =
            BigInt(mantissa);

        if (shift != 0)
            result <<= shift;

        if (representation.sign)
            result = -result;

        return result;
    }


    /*
     * Independent exact oracle for geo3-d's public convention:
     *
     *     sign(det(b-a, c-a, d-a))
     */
    private int oracleOrientation(
        P a,
        P b,
        P c,
        P d
    )
        @safe
    {
        assert(a.isFinite);
        assert(b.isFinite);
        assert(c.isFinite);
        assert(d.isFinite);


        const BigInt ax =
            oracleBinary64Integer(a.x);

        const BigInt ay =
            oracleBinary64Integer(a.y);

        const BigInt az =
            oracleBinary64Integer(a.z);


        const BigInt bx =
            oracleBinary64Integer(b.x);

        const BigInt by =
            oracleBinary64Integer(b.y);

        const BigInt bz =
            oracleBinary64Integer(b.z);


        const BigInt cx =
            oracleBinary64Integer(c.x);

        const BigInt cy =
            oracleBinary64Integer(c.y);

        const BigInt cz =
            oracleBinary64Integer(c.z);


        const BigInt dx =
            oracleBinary64Integer(d.x);

        const BigInt dy =
            oracleBinary64Integer(d.y);

        const BigInt dz =
            oracleBinary64Integer(d.z);


        const BigInt ux =
            bx - ax;

        const BigInt uy =
            by - ay;

        const BigInt uz =
            bz - az;


        const BigInt vx =
            cx - ax;

        const BigInt vy =
            cy - ay;

        const BigInt vz =
            cz - az;


        const BigInt wx =
            dx - ax;

        const BigInt wy =
            dy - ay;

        const BigInt wz =
            dz - az;


        const BigInt determinant =
              ux * vy * wz
            + uy * vz * wx
            + uz * vx * wy
            - uz * vy * wx
            - uy * vx * wz
            - ux * vz * wy;


        if (determinant > 0)
            return 1;

        if (determinant < 0)
            return -1;

        return 0;
    }


    private struct FilterCounts
    {
        size_t negative;
        size_t coplanar;
        size_t positive;
        size_t uncertain;
    }


    private void checkFilter(
        P a,
        P b,
        P c,
        P d,
        ref FilterCounts counts
    )
        @safe
    {
        const int exact =
            oracleOrientation(
                a,
                b,
                c,
                d
            );

        const OrientationFilterResult filtered =
            orientationFilter(
                a,
                b,
                c,
                d
            );


        assert(
            orientationFilterFinite(
                a,
                b,
                c,
                d
            ) ==
            filtered
        );


        final switch (filtered)
        {
            case OrientationFilterResult.negative:
                ++counts.negative;
                break;

            case OrientationFilterResult.coplanar:
                ++counts.coplanar;
                break;

            case OrientationFilterResult.positive:
                ++counts.positive;
                break;

            case OrientationFilterResult.uncertain:
                ++counts.uncertain;
                return;
        }


        /*
         * Core certification invariant:
         *
         * every non-uncertain result must equal the independent exact oracle.
         */
        assert(
            cast(int) filtered ==
            exact
        );
    }


    private ulong nextRandom(ref ulong state)
        pure nothrow @safe @nogc
    {
        assert(state != 0);

        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;

        return state;
    }


    private double randomFiniteDouble(ref ulong state)
        pure nothrow @safe @nogc
    {
        const ulong bits =
            nextRandom(state);

        DoubleRep representation;

        representation.value =
            0.0;

        representation.fraction =
            bits &
            (
                (1UL << 52) -
                1
            );

        ushort exponent =
            cast(ushort)(
                (bits >> 52) &
                0x7ffUL
            );

        if (exponent == 0x7ff)
            exponent = 0x7fe;

        representation.exponent =
            exponent;

        representation.sign =
            (
                bits &
                (1UL << 63)
            ) != 0;

        return representation.value;
    }


    private double randomModerateDouble(ref ulong state)
        pure nothrow @safe @nogc
    {
        return
            cast(double)
            cast(int)
            nextRandom(state);
    }


    private P randomFinitePoint(ref ulong state)
        pure nothrow @safe @nogc
    {
        return P(
            randomFiniteDouble(state),
            randomFiniteDouble(state),
            randomFiniteDouble(state)
        );
    }


    private P randomModeratePoint(ref ulong state)
        pure nothrow @safe @nogc
    {
        return P(
            randomModerateDouble(state),
            randomModerateDouble(state),
            randomModerateDouble(state)
        );
    }
}


@safe unittest
{
    alias P =
        Point3!double;


    const P a =
        P(
            0.0,
            0.0,
            0.0
        );

    const P b =
        P(
            1.0,
            0.0,
            0.0
        );

    const P c =
        P(
            0.0,
            1.0,
            0.0
        );

    const P d =
        P(
            0.0,
            0.0,
            1.0
        );


    /*
     * Canonical geo3-d positive convention.
     */
    assert(
        oracleOrientation(
            a,
            b,
            c,
            d
        ) == 1
    );

    assert(
        orientationFilter(
            a,
            b,
            c,
            d
        ) ==
        OrientationFilterResult.positive
    );


    /*
     * One transposition reverses the certified sign.
     */
    assert(
        orientationFilter(
            a,
            c,
            b,
            d
        ) ==
        OrientationFilterResult.negative
    );


    /*
     * Structurally exact coplanarity.
     */
    const P coplanarD =
        P(
            1.0,
            1.0,
            0.0
        );

    assert(
        oracleOrientation(
            a,
            b,
            c,
            coplanarD
        ) == 0
    );

    assert(
        orientationFilter(
            a,
            b,
            c,
            coplanarD
        ) ==
        OrientationFilterResult.coplanar
    );


    /*
     * Smallest positive binary64 perturbation away from the plane.
     *
     * The filter deliberately refuses subnormal arithmetic.
     */
    enum double minSubnormal =
        0x0.0000000000001p-1022;

    const P nearD =
        P(
            0.25,
            0.25,
            minSubnormal
        );

    assert(
        oracleOrientation(
            a,
            b,
            c,
            nearD
        ) == 1
    );

    assert(
        orientationFilter(
            a,
            b,
            c,
            nearD
        ) ==
        OrientationFilterResult.uncertain
    );


    /*
     * Multiplication underflow with normal operands.
     *
     * Both bdx and cdy are the smallest positive normal binary64 value.
     * Their exact product is non-zero but far below the representable
     * binary64 range and therefore rounds to 0.0.
     *
     * This exercises productUnderflowed(...) directly rather than reaching
     * uncertain through a subnormal coordinate difference.
     */
    const P underflowA =
        P(
            0.0,
            0.0,
            1.0
        );

    const P underflowB =
        P(
            double.min_normal,
            0.0,
            0.0
        );

    const P underflowC =
        P(
            0.0,
            double.min_normal,
            0.0
        );

    const P underflowD =
        P(
            0.0,
            0.0,
            0.0
        );

    assert(
        oracleOrientation(
            underflowA,
            underflowB,
            underflowC,
            underflowD
        ) == -1
    );

    assert(
        orientationFilter(
            underflowA,
            underflowB,
            underflowC,
            underflowD
        ) ==
        OrientationFilterResult.uncertain
    );


    /*
     * Cancellation-heavy near-coplanarity with entirely normal finite
     * intermediates.
     *
     * For the public determinant:
     *
     *     det(b-a, c-a, d-a)
     *
     * the exact value is one binary64 ULP above the coplanar plane:
     *
     *     nextUp(2.0) - 2.0
     *
     * The sign is therefore positive, but the first-stage error bound is
     * larger than the determinant magnitude. The filter must conservatively
     * return uncertain and leave the exact decision to the later fallback.
     */
    enum double aboveTwo =
        0x1.0000000000001p+1;

    const P cancellationA =
        P(
            0.0,
            0.0,
            0.0
        );

    const P cancellationB =
        P(
            1.0,
            0.0,
            1.0
        );

    const P cancellationC =
        P(
            0.0,
            1.0,
            1.0
        );

    const P cancellationD =
        P(
            1.0,
            1.0,
            aboveTwo
        );

    assert(
        oracleOrientation(
            cancellationA,
            cancellationB,
            cancellationC,
            cancellationD
        ) == 1
    );

    assert(
        orientationFilter(
            cancellationA,
            cancellationB,
            cancellationC,
            cancellationD
        ) ==
        OrientationFilterResult.uncertain
    );


    /*
     * Product overflow is never certified.
     */
    const P hugeB =
        P(
            double.max,
            0.0,
            0.0
        );

    const P hugeC =
        P(
            0.0,
            double.max,
            0.0
        );

    const P hugeD =
        P(
            0.0,
            0.0,
            double.max
        );

    assert(
        oracleOrientation(
            a,
            hugeB,
            hugeC,
            hugeD
        ) == 1
    );

    assert(
        orientationFilter(
            a,
            hugeB,
            hugeC,
            hugeD
        ) ==
        OrientationFilterResult.uncertain
    );


    /*
     * Coordinate subtraction itself may overflow although every input
     * coordinate is finite.
     */
    const P overflowA =
        P(
            -double.max,
            0.0,
            0.0
        );

    const P overflowB =
        P(
            0.0,
            1.0,
            0.0
        );

    const P overflowC =
        P(
            0.0,
            0.0,
            1.0
        );

    const P overflowD =
        P(
            double.max,
            0.0,
            0.0
        );

    assert(
        oracleOrientation(
            overflowA,
            overflowB,
            overflowC,
            overflowD
        ) == 1
    );

    assert(
        orientationFilter(
            overflowA,
            overflowB,
            overflowC,
            overflowD
        ) ==
        OrientationFilterResult.uncertain
    );


    /*
     * Non-finite inputs are never certified.
     */
    assert(
        orientationFilter(
            P(
                double.nan,
                0.0,
                0.0
            ),
            b,
            c,
            d
        ) ==
        OrientationFilterResult.uncertain
    );

    assert(
        orientationFilter(
            P(
                double.infinity,
                0.0,
                0.0
            ),
            b,
            c,
            d
        ) ==
        OrientationFilterResult.uncertain
    );
}


@safe unittest
{
    alias P =
        Point3!double;

    ulong state =
        0xbb67_ae85_84ca_a73bUL;

    FilterCounts moderateCounts;
    FilterCounts fullRangeCounts;
    FilterCounts coplanarCounts;

    enum moderateCases =
        10_000;

    enum fullRangeCases =
        2_000;

    enum coplanarCases =
        2_000;


    /*
     * Ordinary-scale cases should be certified by the fast path.
     */
    foreach (_; 0 .. moderateCases)
    {
        checkFilter(
            randomModeratePoint(state),
            randomModeratePoint(state),
            randomModeratePoint(state),
            randomModeratePoint(state),
            moderateCounts
        );
    }

    assert(
        moderateCounts.uncertain ==
        0
    );


    /*
     * Arbitrary finite binary64 patterns exercise exponent imbalance,
     * overflow, underflow, and conservative fallback.
     *
     * Every result that is certified is checked against the exact oracle.
     */
    foreach (_; 0 .. fullRangeCases)
    {
        checkFilter(
            randomFinitePoint(state),
            randomFinitePoint(state),
            randomFinitePoint(state),
            randomFinitePoint(state),
            fullRangeCounts
        );
    }

    assert(
        fullRangeCounts.uncertain >
        0
    );


    /*
     * Exact z=0 geometry is structurally coplanar.
     */
    foreach (_; 0 .. coplanarCases)
    {
        const P ca =
            P(
                randomModerateDouble(state),
                randomModerateDouble(state),
                0.0
            );

        const P cb =
            P(
                randomModerateDouble(state),
                randomModerateDouble(state),
                0.0
            );

        const P cc =
            P(
                randomModerateDouble(state),
                randomModerateDouble(state),
                0.0
            );

        const P cd =
            P(
                randomModerateDouble(state),
                randomModerateDouble(state),
                0.0
            );

        checkFilter(
            ca,
            cb,
            cc,
            cd,
            coplanarCounts
        );
    }

    assert(
        coplanarCounts.coplanar ==
        coplanarCases
    );

    assert(
        coplanarCounts.uncertain ==
        0
    );
}
