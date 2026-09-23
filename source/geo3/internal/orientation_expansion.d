/**
 * Exact binary64 expansion fallback for three-dimensional orientation.
 *
 * INTERNAL IMPLEMENTATION MODULE.
 *
 * This backend covers a deliberately conservative ordinary exponent range.
 * Inputs outside its validated working range return false and belong to the
 * later full-range dyadic fallback.
 *
 * No tolerance or epsilon participates in the result.
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
module geo3.internal.orientation_expansion;

import geo3.internal.binary64_rounding :
    roundedSub;

import geo3.internal.expansion :
    ExpansionBuffer,
    TwoComponent,
    expansionSign,
    fastExpansionSumZeroElim,
    negateExpansion,
    scaleExpansionZeroElim,
    twoDiff;

import geo3.point :
    Point3;

import std.math.traits :
    isFinite;


/*
 * Conservative validated expansion working range.
 *
 * These limits keep splitting, products, triple products, and their
 * roundoff components well inside the normal finite binary64 range.
 *
 * This is not the complete finite binary64 domain. Slice 4 provides the
 * exact full-range dyadic backend.
 */
private enum double minWorkingMagnitude =
    0x1p-200;

private enum double maxWorkingMagnitude =
    0x1p+200;


/*
 * Absolute magnitude for a known finite value.
 */
private double magnitude(double value)
    pure nothrow @safe @nogc
{
    return value < 0.0
        ? -value
        : value;
}


/*
 * True when one expansion component belongs to the validated working range.
 *
 * Exact zero is always supported.
 */
private bool supportedComponent(double value)
    pure nothrow @safe @nogc
{
    if (!isFinite(value))
        return false;

    if (value == 0.0)
        return true;

    const double absolute =
        magnitude(value);

    return
        absolute >= minWorkingMagnitude &&
        absolute <= maxWorkingMagnitude;
}


/*
 * Builds the exact expansion of lhs-rhs.
 *
 * Components are stored least-significant first.
 *
 * false means that the difference lies outside the conservative expansion
 * working range and must be handled by the full-range dyadic backend.
 */
private bool buildDifference(
    double lhs,
    double rhs,
    ref ExpansionBuffer!2 result
)
    pure nothrow @safe @nogc
{
    const double rounded =
        roundedSub(
            lhs,
            rhs
        );

    if (!isFinite(rounded))
        return false;


    const TwoComponent difference =
        twoDiff(
            lhs,
            rhs
        );


    if (
        !supportedComponent(difference.low) ||
        !supportedComponent(difference.high)
    )
        return false;


    result.clear();

    if (difference.low != 0.0)
    {
        result.append(
            difference.low
        );
    }

    if (difference.high != 0.0)
    {
        result.append(
            difference.high
        );
    }

    if (result.empty)
        result.append(0.0);


    return true;
}


/*
 * Exact multiplication by an at-most-two-component coordinate difference.
 */
private void multiplyByDifference(
    size_t SourceCapacity,
    size_t ResultCapacity
)(
    ref const ExpansionBuffer!SourceCapacity source,
    ref const ExpansionBuffer!2 difference,
    ref ExpansionBuffer!ResultCapacity result
)
    pure nothrow @safe @nogc
if (ResultCapacity >= 4 * SourceCapacity)
{
    assert(!source.empty);
    assert(!difference.empty);
    assert(difference.length <= 2);


    ExpansionBuffer!(2 * SourceCapacity) first;
    ExpansionBuffer!(2 * SourceCapacity) second;


    scaleExpansionZeroElim(
        source,
        difference[0],
        first
    );


    if (difference.length == 2)
    {
        scaleExpansionZeroElim(
            source,
            difference[1],
            second
        );
    }


    fastExpansionSumZeroElim(
        first,
        second,
        result
    );
}


/*
 * Exact product of three at-most-two-component coordinate differences.
 *
 * Capacity derivation:
 *
 *     2 x 2 -> <= 8 components
 *     8 x 2 -> <= 32 components
 */
private void tripleProduct(
    ref const ExpansionBuffer!2 x,
    ref const ExpansionBuffer!2 y,
    ref const ExpansionBuffer!2 z,
    ref ExpansionBuffer!32 result
)
    pure nothrow @safe @nogc
{
    ExpansionBuffer!8 xy;

    multiplyByDifference(
        x,
        y,
        xy
    );

    multiplyByDifference(
        xy,
        z,
        result
    );
}


/**
 * Attempts exact 3D orientation through binary64 expansion arithmetic.
 *
 * The represented determinant is exactly:
 *
 *     det(b-a, c-a, d-a)
 *
 * On success:
 *
 *     sign = -1, 0, or 1
 *
 * and the returned sign is mathematically exact.
 *
 * false does not mean degeneracy. It means only that the finite input lies
 * outside this conservative expansion backend's validated working range and
 * must be routed to the later full-range dyadic backend.
 *
 * Non-finite inputs also return false. The eventual public predicate retains
 * its separate finite-coordinate precondition.
 */
bool tryOrientationExactExpansion(
    Point3!double a,
    Point3!double b,
    Point3!double c,
    Point3!double d,
    out int sign
)
    pure nothrow @safe @nogc
{
    if (
        !a.isFinite ||
        !b.isFinite ||
        !c.isFinite ||
        !d.isFinite
    )
        return false;


    ExpansionBuffer!2 ux;
    ExpansionBuffer!2 uy;
    ExpansionBuffer!2 uz;

    ExpansionBuffer!2 vx;
    ExpansionBuffer!2 vy;
    ExpansionBuffer!2 vz;

    ExpansionBuffer!2 wx;
    ExpansionBuffer!2 wy;
    ExpansionBuffer!2 wz;


    if (
        !buildDifference(
            b.x,
            a.x,
            ux
        ) ||
        !buildDifference(
            b.y,
            a.y,
            uy
        ) ||
        !buildDifference(
            b.z,
            a.z,
            uz
        ) ||
        !buildDifference(
            c.x,
            a.x,
            vx
        ) ||
        !buildDifference(
            c.y,
            a.y,
            vy
        ) ||
        !buildDifference(
            c.z,
            a.z,
            vz
        ) ||
        !buildDifference(
            d.x,
            a.x,
            wx
        ) ||
        !buildDifference(
            d.y,
            a.y,
            wy
        ) ||
        !buildDifference(
            d.z,
            a.z,
            wz
        )
    )
        return false;


    ExpansionBuffer!32 t0;
    ExpansionBuffer!32 t1;
    ExpansionBuffer!32 t2;

    ExpansionBuffer!32 raw3;
    ExpansionBuffer!32 raw4;
    ExpansionBuffer!32 raw5;

    ExpansionBuffer!32 t3;
    ExpansionBuffer!32 t4;
    ExpansionBuffer!32 t5;


    /*
     * det(u, v, w)
     *
     *   + ux * vy * wz
     *   + uy * vz * wx
     *   + uz * vx * wy
     *   - uz * vy * wx
     *   - uy * vx * wz
     *   - ux * vz * wy
     */
    tripleProduct(
        ux,
        vy,
        wz,
        t0
    );

    tripleProduct(
        uy,
        vz,
        wx,
        t1
    );

    tripleProduct(
        uz,
        vx,
        wy,
        t2
    );

    tripleProduct(
        uz,
        vy,
        wx,
        raw3
    );

    tripleProduct(
        uy,
        vx,
        wz,
        raw4
    );

    tripleProduct(
        ux,
        vz,
        wy,
        raw5
    );


    negateExpansion(
        raw3,
        t3
    );

    negateExpansion(
        raw4,
        t4
    );

    negateExpansion(
        raw5,
        t5
    );


    /*
     * Balanced exact sum tree.
     */
    ExpansionBuffer!64 s01;
    ExpansionBuffer!64 s23;
    ExpansionBuffer!64 s45;


    fastExpansionSumZeroElim(
        t0,
        t1,
        s01
    );

    fastExpansionSumZeroElim(
        t2,
        t3,
        s23
    );

    fastExpansionSumZeroElim(
        t4,
        t5,
        s45
    );


    ExpansionBuffer!128 s0123;

    fastExpansionSumZeroElim(
        s01,
        s23,
        s0123
    );


    ExpansionBuffer!192 determinant;

    fastExpansionSumZeroElim(
        s0123,
        s45,
        determinant
    );


    sign =
        expansionSign(
            determinant
        );


    return true;
}


version(unittest)
{
    import geo3.internal.orientation_filter :
        OrientationFilterResult,
        orientationFilter;

    import std.bigint :
        BigInt;

    import std.bitmanip :
        DoubleRep;


    private alias P =
        Point3!double;


    /*
     * Exact finite binary64 decoding as an integer multiple of 2^-1074.
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

            shift =
                0;
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
     * Independent exact oracle for:
     *
     *     det(b-a, c-a, d-a)
     */
    private int oracleOrientation(
        P a,
        P b,
        P c,
        P d
    )
        @safe
    {
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


    private void requireExpansion(
        P a,
        P b,
        P c,
        P d
    )
        @safe
    {
        int sign;

        assert(
            tryOrientationExactExpansion(
                a,
                b,
                c,
                d,
                sign
            )
        );


        const int expected =
            oracleOrientation(
                a,
                b,
                c,
                d
            );


        assert(
            sign ==
            expected
        );


        int reversed;

        assert(
            tryOrientationExactExpansion(
                a,
                c,
                b,
                d,
                reversed
            )
        );

        assert(
            reversed ==
            -sign
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


    private double randomModerate(ref ulong state)
        pure nothrow @safe @nogc
    {
        return
            cast(double)
            cast(int)
            nextRandom(state);
    }


    private P randomModeratePoint(ref ulong state)
        pure nothrow @safe @nogc
    {
        return P(
            randomModerate(state),
            randomModerate(state),
            randomModerate(state)
        );
    }
}


@safe unittest
{
    alias P =
        Point3!double;


    static assert(
        supportedComponent(
            0.0
        )
    );

    static assert(
        supportedComponent(
            0x1p-200
        )
    );

    static assert(
        supportedComponent(
            0x1p+200
        )
    );

    static assert(
        !supportedComponent(
            0x1p-201
        )
    );

    static assert(
        !supportedComponent(
            0x1p+201
        )
    );


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
            1.0
        );

    const P c =
        P(
            0.0,
            1.0,
            1.0
        );


    /*
     * Cancellation-heavy exact coplanarity.
     */
    requireExpansion(
        a,
        b,
        c,
        P(
            1.0,
            1.0,
            2.0
        )
    );


    /*
     * One ULP above the coplanar plane.
     *
     * This is uncertain in the first-stage filter but exactly positive in
     * the expansion fallback.
     */
    const P above =
        P(
            1.0,
            1.0,
            0x1.0000000000001p+1
        );

    assert(
        orientationFilter(
            a,
            b,
            c,
            above
        ) ==
        OrientationFilterResult.uncertain
    );

    requireExpansion(
        a,
        b,
        c,
        above
    );


    /*
     * One representable binary64 value below 2.0.
     *
     * Exact sign is the opposite side of the same plane.
     */
    const P below =
        P(
            1.0,
            1.0,
            0x1.fffffffffffffp+0
        );

    assert(
        oracleOrientation(
            a,
            b,
            c,
            below
        ) ==
        -1
    );

    requireExpansion(
        a,
        b,
        c,
        below
    );


    /*
     * Exercise a genuine non-zero TwoDiff tail.
     */
    requireExpansion(
        P(
            0x1p-100,
            0.0,
            0.0
        ),
        P(
            1.0,
            0.0,
            0.0
        ),
        P(
            0.0,
            1.0,
            0.0
        ),
        P(
            0.0,
            0.0,
            1.0
        )
    );


    /*
     * Complete degeneracy.
     */
    requireExpansion(
        P(
            7.0,
            -2.0,
            9.0
        ),
        P(
            7.0,
            -2.0,
            9.0
        ),
        P(
            7.0,
            -2.0,
            9.0
        ),
        P(
            7.0,
            -2.0,
            9.0
        )
    );


    /*
     * Even extreme coordinates are harmless when every exact difference
     * is zero.
     */
    requireExpansion(
        P(
            double.max,
            double.max,
            double.max
        ),
        P(
            double.max,
            double.max,
            double.max
        ),
        P(
            double.max,
            double.max,
            double.max
        ),
        P(
            double.max,
            double.max,
            double.max
        )
    );


    int unusedSign;


    /*
     * Finite coordinate subtraction overflow belongs to Slice 4.
     */
    assert(
        !tryOrientationExactExpansion(
            P(
                -double.max,
                0.0,
                0.0
            ),
            P(
                double.max,
                0.0,
                0.0
            ),
            P(
                0.0,
                1.0,
                0.0
            ),
            P(
                0.0,
                0.0,
                1.0
            ),
            unusedSign
        )
    );


    /*
     * Subnormal-scale geometry likewise falls through to the full-range
     * dyadic backend.
     */
    enum double minSubnormal =
        0x0.0000000000001p-1022;

    const P subnormalB =
        P(
            minSubnormal,
            0.0,
            0.0
        );

    const P subnormalC =
        P(
            0.0,
            minSubnormal,
            0.0
        );

    const P subnormalD =
        P(
            0.0,
            0.0,
            minSubnormal
        );

    assert(
        orientationFilter(
            P.init,
            subnormalB,
            subnormalC,
            subnormalD
        ) ==
        OrientationFilterResult.uncertain
    );

    assert(
        !tryOrientationExactExpansion(
            P.init,
            subnormalB,
            subnormalC,
            subnormalD,
            unusedSign
        )
    );


    /*
     * Non-finite input cannot enter the exact expansion backend.
     */
    assert(
        !tryOrientationExactExpansion(
            P(
                double.nan,
                0.0,
                0.0
            ),
            b,
            c,
            P(
                0.0,
                0.0,
                1.0
            ),
            unusedSign
        )
    );
}


@safe unittest
{
    /*
     * Production regression sweep.
     *
     * The research prototype already validated 50,000 deterministic
     * moderate-range quadruples under DMD/LDC debug/release. Keep a smaller
     * deterministic exact-oracle sweep in the normal production test suite.
     */
    ulong state =
        0x3c6e_f372_fe94_f82bUL;

    enum randomCases =
        10_000;


    foreach (_; 0 .. randomCases)
    {
        requireExpansion(
            randomModeratePoint(state),
            randomModeratePoint(state),
            randomModeratePoint(state),
            randomModeratePoint(state)
        );
    }
}
