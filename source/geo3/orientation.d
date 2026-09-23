/**
 * Robust orientation predicates for three-dimensional points.
 *
 * Orientation3 classifies four points by the exact sign of:
 *
 *     det(b-a, c-a, d-a)
 *
 * No tolerance or epsilon participates in the predicate.
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
module geo3.orientation;

import geo3.internal.orientation_dyadic :
    orientationDyadicExact;

import geo3.internal.orientation_expansion :
    tryOrientationExactExpansion;

import geo3.internal.orientation_filter :
    OrientationFilterResult,
    orientationFilter;

import geo3.internal.orientation_integral :
    orientationIntegralSign;

import geo3.point :
    Point3;


/**
 * Exact three-dimensional orientation / coplanarity result.
 *
 * The sign convention is defined by:
 *
 *     det(b-a, c-a, d-a)
 *
 * `Orientation3.init` is deliberately the neutral `coplanar` state.
 */
enum Orientation3 : byte
{
    /// The determinant is exactly zero.
    coplanar = 0,

    /// The determinant is negative.
    negative = -1,

    /// The determinant is positive.
    positive = 1,
}


/*
 * Converts an exact determinant sign into the public enum.
 */
private Orientation3 fromDeterminantSign(int sign)
    pure nothrow @safe @nogc
{
    if (sign > 0)
        return Orientation3.positive;

    if (sign < 0)
        return Orientation3.negative;

    return Orientation3.coplanar;
}


/**
 * Exact orientation predicate for Point3!int.
 *
 * Returns the mathematically exact sign of:
 *
 *     det(b-a, c-a, d-a)
 *
 * over the complete int coordinate domain.
 *
 * Degenerate affine configurations are valid inputs. Repeated points,
 * collinear triples, and all other exact zero-determinant configurations
 * return `Orientation3.coplanar`.
 *
 * No allocation is performed.
 *
 * Complexity:
 *     O(1) time and O(1) auxiliary space.
 */
Orientation3 orientation(
    Point3!int a,
    Point3!int b,
    Point3!int c,
    Point3!int d
)
    pure nothrow @safe @nogc
{
    return fromDeterminantSign(
        orientationIntegralSign(
            a,
            b,
            c,
            d
        )
    );
}


/**
 * Exact orientation predicate for Point3!long.
 *
 * Returns the mathematically exact sign of:
 *
 *     det(b-a, c-a, d-a)
 *
 * over the complete long coordinate domain.
 *
 * Intermediate determinant terms exceed native 64- and 128-bit signed
 * arithmetic at the coordinate extremes. The implementation therefore uses
 * fixed-width exact unsigned arithmetic internally and determines only the
 * final determinant sign.
 *
 * Degenerate affine configurations are valid inputs.
 *
 * No allocation is performed.
 *
 * Complexity:
 *     O(1) time and O(1) auxiliary space.
 */
Orientation3 orientation(
    Point3!long a,
    Point3!long b,
    Point3!long c,
    Point3!long d
)
    pure nothrow @safe @nogc
{
    return fromDeterminantSign(
        orientationIntegralSign(
            a,
            b,
            c,
            d
        )
    );
}

/**
 * Robust exact orientation predicate for Point3!double.
 *
 * Preconditions:
 *
 *     all coordinates are finite.
 *
 * Returns the mathematically exact sign of:
 *
 *     det(b-a, c-a, d-a)
 *
 * over the complete finite IEEE binary64 coordinate domain.
 *
 * The implementation uses a staged backend:
 *
 * 1. a certified binary64 floating-point filter;
 * 2. exact expansion arithmetic in its validated working range;
 * 3. an exact fixed-width dyadic fallback for every remaining finite input.
 *
 * The internal `uncertain` filter state is never exposed through the public
 * API.
 *
 * Degenerate affine configurations are valid inputs and return
 * `Orientation3.coplanar` exactly when the determinant is mathematically zero.
 *
 * No allocation is performed.
 *
 * Complexity:
 *     O(1) time and O(1) auxiliary space.
 */
Orientation3 orientation(
    Point3!double a,
    Point3!double b,
    Point3!double c,
    Point3!double d
)
    pure nothrow @safe @nogc
{
    assert(a.isFinite);
    assert(b.isFinite);
    assert(c.isFinite);
    assert(d.isFinite);


    const OrientationFilterResult filtered =
        orientationFilter(
            a,
            b,
            c,
            d
        );


    final switch (filtered)
    {
        case OrientationFilterResult.negative:
            return Orientation3.negative;

        case OrientationFilterResult.coplanar:
            return Orientation3.coplanar;

        case OrientationFilterResult.positive:
            return Orientation3.positive;

        case OrientationFilterResult.uncertain:
            break;
    }


    int sign;


    if (
        tryOrientationExactExpansion(
            a,
            b,
            c,
            d,
            sign
        )
    )
    {
        return fromDeterminantSign(
            sign
        );
    }


    /*
     * Every finite binary64 input not covered by the expansion backend is
     * exactly representable by the full-range fixed-width dyadic backend.
     */
    return fromDeterminantSign(
        orientationDyadicExact(
            a,
            b,
            c,
            d
        )
    );
}



static assert(
    Orientation3.init ==
    Orientation3.coplanar
);

static assert(
    cast(byte) Orientation3.negative == -1
);

static assert(
    cast(byte) Orientation3.coplanar == 0
);

static assert(
    cast(byte) Orientation3.positive == 1
);


/*
 * binary64 is supported by the complete filter -> expansion -> dyadic
 * pipeline.
 *
 * binary32 remains deliberately deferred to Slice 5, where exact promotion
 * to binary64 will reuse this backend.
 *
 * real remains unsupported until a platform-aware robust backend is designed.
 */
static assert(
    !__traits(
        compiles,
        orientation(
            Point3!float.init,
            Point3!float.init,
            Point3!float.init,
            Point3!float.init
        )
    )
);

static assert(
    __traits(
        compiles,
        orientation(
            Point3!double.init,
            Point3!double.init,
            Point3!double.init,
            Point3!double.init
        )
    )
);

static assert(
    !__traits(
        compiles,
        orientation(
            Point3!real.init,
            Point3!real.init,
            Point3!real.init,
            Point3!real.init
        )
    )
);


version(unittest)
{
    import std.bigint :
        BigInt;

    import std.bitmanip :
        DoubleRep;

    import std.math.traits :
        isFinite;


    /*
     * Independent arbitrary-precision oracle.
     *
     * Production orientation does not depend on BigInt.
     */
    private Orientation3 oracleOrientation(T)(
        Point3!T a,
        Point3!T b,
        Point3!T c,
        Point3!T d
    )
        @safe
    if (
        is(T == int) ||
        is(T == long)
    )
    {
        const BigInt ux =
            BigInt(b.x) -
            BigInt(a.x);

        const BigInt uy =
            BigInt(b.y) -
            BigInt(a.y);

        const BigInt uz =
            BigInt(b.z) -
            BigInt(a.z);


        const BigInt vx =
            BigInt(c.x) -
            BigInt(a.x);

        const BigInt vy =
            BigInt(c.y) -
            BigInt(a.y);

        const BigInt vz =
            BigInt(c.z) -
            BigInt(a.z);


        const BigInt wx =
            BigInt(d.x) -
            BigInt(a.x);

        const BigInt wy =
            BigInt(d.y) -
            BigInt(a.y);

        const BigInt wz =
            BigInt(d.z) -
            BigInt(a.z);


        const BigInt determinant =
              ux * vy * wz
            + uy * vz * wx
            + uz * vx * wy
            - uz * vy * wx
            - uy * vx * wz
            - ux * vz * wy;


        if (determinant > 0)
            return Orientation3.positive;

        if (determinant < 0)
            return Orientation3.negative;

        return Orientation3.coplanar;
    }

    /*
     * Independent exact finite binary64 decoding:
     *
     *     value = integer * 2^-1074
     *
     * Production robust-predicate arithmetic is not reused here.
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
     * Independent arbitrary-precision oracle for Point3!double.
     */
    private Orientation3 oracleDoubleOrientation(
        Point3!double a,
        Point3!double b,
        Point3!double c,
        Point3!double d
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
            return Orientation3.positive;

        if (determinant < 0)
            return Orientation3.negative;

        return Orientation3.coplanar;
    }



    private Orientation3 opposite(
        Orientation3 value
    )
        pure nothrow @safe @nogc
    {
        final switch (value)
        {
            case Orientation3.negative:
                return Orientation3.positive;

            case Orientation3.coplanar:
                return Orientation3.coplanar;

            case Orientation3.positive:
                return Orientation3.negative;
        }
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


    private Point3!int randomIntPoint(ref ulong state)
        pure nothrow @safe @nogc
    {
        return Point3!int(
            cast(int) nextRandom(state),
            cast(int) nextRandom(state),
            cast(int) nextRandom(state)
        );
    }


    private Point3!long randomLongPoint(ref ulong state)
        pure nothrow @safe @nogc
    {
        return Point3!long(
            cast(long) nextRandom(state),
            cast(long) nextRandom(state),
            cast(long) nextRandom(state)
        );
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
                (
                    bits >>
                    52
                ) &
                0x7ffUL
            );


        /*
         * Raw exponent 2047 represents infinity or NaN.
         * Remap it to the largest finite exponent.
         */
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


    private Point3!double randomDoublePoint(ref ulong state)
        pure nothrow @safe @nogc
    {
        return Point3!double(
            randomFiniteDouble(state),
            randomFiniteDouble(state),
            randomFiniteDouble(state)
        );
    }


    private void verifyDoubleAgainstOracle(
        Point3!double a,
        Point3!double b,
        Point3!double c,
        Point3!double d
    )
        @safe
    {
        const Orientation3 expected =
            oracleDoubleOrientation(
                a,
                b,
                c,
                d
            );

        const Orientation3 actual =
            orientation(
                a,
                b,
                c,
                d
            );


        assert(
            actual ==
            expected
        );


        assert(
            orientation(
                a,
                c,
                b,
                d
            ) ==
            opposite(actual)
        );


        /*
         * A four-cycle is odd.
         */
        assert(
            orientation(
                b,
                c,
                d,
                a
            ) ==
            opposite(actual)
        );


        /*
         * Two transpositions preserve sign.
         */
        assert(
            orientation(
                b,
                a,
                d,
                c
            ) ==
            actual
        );
    }



    private void verifyAgainstOracle(T)(
        Point3!T a,
        Point3!T b,
        Point3!T c,
        Point3!T d
    )
        @safe
    if (
        is(T == int) ||
        is(T == long)
    )
    {
        const Orientation3 expected =
            oracleOrientation(
                a,
                b,
                c,
                d
            );

        const Orientation3 actual =
            orientation(
                a,
                b,
                c,
                d
            );

        assert(actual == expected);


        /*
         * Any single transposition reverses sign.
         */
        assert(
            orientation(
                a,
                c,
                b,
                d
            ) ==
            opposite(actual)
        );


        /*
         * The four-cycle (a,b,c,d) -> (b,c,d,a) is odd.
         */
        assert(
            orientation(
                b,
                c,
                d,
                a
            ) ==
            opposite(actual)
        );


        /*
         * Two transpositions preserve sign.
         */
        assert(
            orientation(
                b,
                a,
                d,
                c
            ) ==
            actual
        );
    }
}


@safe unittest
{
    alias PI =
        Point3!int;

    alias PL =
        Point3!long;


    /*
     * Canonical public sign convention.
     */
    assert(
        orientation(
            PI(0, 0, 0),
            PI(1, 0, 0),
            PI(0, 1, 0),
            PI(0, 0, 1)
        ) ==
        Orientation3.positive
    );


    /*
     * One transposition reverses the canonical sign.
     */
    assert(
        orientation(
            PI(0, 0, 0),
            PI(0, 1, 0),
            PI(1, 0, 0),
            PI(0, 0, 1)
        ) ==
        Orientation3.negative
    );


    /*
     * Exact coplanarity.
     */
    assert(
        orientation(
            PI(0, 0, 0),
            PI(1, 0, 0),
            PI(0, 1, 0),
            PI(1, 1, 0)
        ) ==
        Orientation3.coplanar
    );


    /*
     * Repeated-point degeneracy.
     */
    assert(
        orientation(
            PI(7, -2, 9),
            PI(7, -2, 9),
            PI(3, 5, 11),
            PI(-8, 4, 6)
        ) ==
        Orientation3.coplanar
    );


    /*
     * Complete degeneracy.
     */
    assert(
        orientation(
            PI(4, 5, 6),
            PI(4, 5, 6),
            PI(4, 5, 6),
            PI(4, 5, 6)
        ) ==
        Orientation3.coplanar
    );


    /*
     * Collinear first three points imply zero determinant for every d.
     */
    assert(
        orientation(
            PI(0, 0, 0),
            PI(1, 2, 3),
            PI(2, 4, 6),
            PI(9, -7, 5)
        ) ==
        Orientation3.coplanar
    );


    /*
     * Bucket-width stress case for int.
     *
     * Relative to a, the component-sign pattern is:
     *
     *     u = (-, -, -)
     *     v = (-, -, +)
     *     w = (-, +, -)
     *
     * The six Leibniz terms therefore split 5 positive / 1 negative.
     *
     * With the complete X span and half-range Y/Z spans used here, the
     * positive exact bucket requires 97 bits. A three-limb (96-bit) bucket
     * would overflow; the validated fourth limb is genuinely required.
     */
    verifyAgainstOracle(
        PI(
            int.max,
            0,
            0
        ),
        PI(
            int.min,
            int.min,
            int.min
        ),
        PI(
            int.min,
            int.min,
            int.max
        ),
        PI(
            int.min,
            int.max,
            int.min
        )
    );


    /*
     * Bucket-width stress case for long.
     *
     * This is the 64-bit integral analogue of the int case above.
     *
     * Five determinant terms accumulate in the positive bucket. Its exact
     * magnitude requires 193 bits, so six uint limbs (192 bits) would be
     * insufficient. The validated seventh limb is genuinely required.
     */
    verifyAgainstOracle(
        PL(
            long.max,
            0,
            0
        ),
        PL(
            long.min,
            long.min,
            long.min
        ),
        PL(
            long.min,
            long.min,
            long.max
        ),
        PL(
            long.min,
            long.max,
            long.min
        )
    );


    /*
     * Complete int coordinate span.
     */
    verifyAgainstOracle(
        PI(
            int.min,
            int.min,
            int.min
        ),
        PI(
            int.max,
            int.min,
            int.min
        ),
        PI(
            int.min,
            int.max,
            int.min
        ),
        PI(
            int.min,
            int.min,
            int.max
        )
    );


    /*
     * Complete long coordinate span.
     *
     * The determinant magnitude is approximately (2^64)^3.
     */
    verifyAgainstOracle(
        PL(
            long.min,
            long.min,
            long.min
        ),
        PL(
            long.max,
            long.min,
            long.min
        ),
        PL(
            long.min,
            long.max,
            long.min
        ),
        PL(
            long.min,
            long.min,
            long.max
        )
    );


    /*
     * Exact coplanarity across complete X/Y spans.
     */
    verifyAgainstOracle(
        PL(
            long.min,
            long.min,
            long.min
        ),
        PL(
            long.max,
            long.min,
            long.min
        ),
        PL(
            long.min,
            long.max,
            long.min
        ),
        PL(
            long.max,
            long.max,
            long.min
        )
    );


    /*
     * Large cancellation-heavy coplanarity with no zero direction
     * components.
     */
    enum long m =
        (1L << 61) - 1;

    verifyAgainstOracle(
        PL(0, 0, 0),
        PL(
            m,
            2 * m,
            3 * m
        ),
        PL(
            2 * m,
            -m,
            m
        ),
        PL(
            3 * m,
            m,
            4 * m
        )
    );


    /*
     * Deterministic exact-oracle sweeps across full integral bit patterns.
     */
    ulong state =
        0x6a09_e667_f3bc_c909UL;

    enum randomCases =
        10_000;

    foreach (_; 0 .. randomCases)
    {
        verifyAgainstOracle(
            randomIntPoint(state),
            randomIntPoint(state),
            randomIntPoint(state),
            randomIntPoint(state)
        );
    }

    foreach (_; 0 .. randomCases)
    {
        verifyAgainstOracle(
            randomLongPoint(state),
            randomLongPoint(state),
            randomLongPoint(state),
            randomLongPoint(state)
        );
    }
}


@safe unittest
{
    alias P =
        Point3!double;


    const P zero =
        P(
            0.0,
            0.0,
            0.0
        );


    /*
     * Canonical public binary64 sign convention.
     */
    verifyDoubleAgainstOracle(
        zero,
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

    assert(
        orientation(
            zero,
            P(1.0, 0.0, 0.0),
            P(0.0, 1.0, 0.0),
            P(0.0, 0.0, 1.0)
        ) ==
        Orientation3.positive
    );


    /*
     * Near-coplanar case routed beyond the first-stage filter.
     */
    verifyDoubleAgainstOracle(
        zero,
        P(
            1.0,
            0.0,
            1.0
        ),
        P(
            0.0,
            1.0,
            1.0
        ),
        P(
            1.0,
            1.0,
            0x1.0000000000001p+1
        )
    );


    /*
     * Smallest positive binary64 subnormal in all three basis directions.
     *
     * This necessarily exercises the full-range dyadic fallback.
     */
    enum double minSubnormal =
        0x0.0000000000001p-1022;

    verifyDoubleAgainstOracle(
        zero,
        P(
            minSubnormal,
            0.0,
            0.0
        ),
        P(
            0.0,
            minSubnormal,
            0.0
        ),
        P(
            0.0,
            0.0,
            minSubnormal
        )
    );


    /*
     * Coordinate subtraction itself exceeds binary64.
     */
    verifyDoubleAgainstOracle(
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
            -double.max,
            1.0,
            0.0
        ),
        P(
            -double.max,
            0.0,
            1.0
        )
    );


    /*
     * Full-span exact coplanarity.
     */
    verifyDoubleAgainstOracle(
        P(
            -double.max,
            -double.max,
            0.0
        ),
        P(
            double.max,
            -double.max,
            0.0
        ),
        P(
            -double.max,
            double.max,
            0.0
        ),
        P(
            double.max,
            double.max,
            0.0
        )
    );


    /*
     * Five-versus-one determinant-term sign pattern at extreme binary64
     * magnitudes. This exercises repeated exact accumulation into one
     * 198-limb determinant bucket.
     */
    verifyDoubleAgainstOracle(
        P(
            double.max,
            0.0,
            0.0
        ),
        P(
            -double.max,
            -double.max,
            -double.max
        ),
        P(
            -double.max,
            -double.max,
            double.max
        ),
        P(
            -double.max,
            double.max,
            -double.max
        )
    );


    /*
     * Deterministic complete-finite-binary64 public-pipeline sweep.
     *
     * The independent BigInt oracle decides the exact result while the
     * production function chooses among filter, expansion, and dyadic paths.
     */
    ulong state =
        0x243f_6a88_85a3_08d3UL;

    enum randomCases =
        512;


    foreach (_; 0 .. randomCases)
    {
        verifyDoubleAgainstOracle(
            randomDoublePoint(state),
            randomDoublePoint(state),
            randomDoublePoint(state),
            randomDoublePoint(state)
        );
    }
}
