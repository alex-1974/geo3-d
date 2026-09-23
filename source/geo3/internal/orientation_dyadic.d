/**
 * Full-range exact binary64 backend for robust 3D orientation.
 *
 * INTERNAL IMPLEMENTATION MODULE.
 *
 * Every finite binary64 coordinate is decoded exactly into a fixed-width
 * dyadic representation. After decoding, no floating-point arithmetic is
 * performed.
 *
 * Numerical widths:
 *
 *     coordinate / difference:  66 limbs
 *     pair product:            132 limbs
 *     triple product / bucket: 198 limbs
 *
 * 198 * 32 = 6336 bits.
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
module geo3.internal.orientation_dyadic;

import geo3.internal.dyadic :
    SignedDyadicDifference,
    decodeBinary64Coordinate,
    dyadicCoordinateLimbs,
    dyadicProductLimbs,
    multiplyDyadicDifferences,
    subtractDyadicCoordinates;

import geo3.internal.fixed_uint :
    UIntFixed,
    addAssign,
    compareUnsigned,
    multiplyUnsigned;

import geo3.point :
    Point3;

import std.math.traits :
    isFinite;


/**
 * Number of 32-bit limbs required by every exact triple product and complete
 * positive/negative Orientation3 determinant bucket.
 */
enum size_t orientationTripleLimbs =
    3 * dyadicCoordinateLimbs;


alias OrientationTripleMagnitude =
    UIntFixed!orientationTripleLimbs;


static assert(
    dyadicCoordinateLimbs ==
    66
);

static assert(
    dyadicProductLimbs ==
    132
);

static assert(
    orientationTripleLimbs ==
    198
);

static assert(
    orientationTripleLimbs * 32 ==
    6336
);


/*
 * Accumulates one exact signed Leibniz triple product into its unsigned
 * determinant-sign bucket.
 */
private void accumulateTriple(
    ref OrientationTripleMagnitude positive,
    ref OrientationTripleMagnitude negative,
    int coefficient,
    ref const SignedDyadicDifference x,
    ref const SignedDyadicDifference y,
    ref const SignedDyadicDifference z
)
    pure nothrow @safe @nogc
{
    assert(
        coefficient == -1 ||
        coefficient == 1
    );


    if (
        x.sign == 0 ||
        y.sign == 0 ||
        z.sign == 0
    )
        return;


    const auto xy =
        multiplyDyadicDifferences(
            x,
            y
        );

    assert(
        xy.sign != 0
    );


    const OrientationTripleMagnitude xyz =
        multiplyUnsigned(
            xy.magnitude,
            z.magnitude
        );


    const int sign =
        coefficient *
        xy.sign *
        z.sign;


    const bool success =
        sign > 0
            ? addAssign(
                positive,
                xyz
            )
            : addAssign(
                negative,
                xyz
            );


    /*
     * Maximum coordinate-difference magnitude is below 2^2099.
     *
     * Therefore each triple product is below:
     *
     *     2^6297
     *
     * At most five non-zero Leibniz terms can share one sign, so one
     * complete sign bucket is below:
     *
     *     5 * 2^6297 < 2^6300
     *
     * 198 limbs provide 6336 bits.
     */
    assert(success);
}


/**
 * Exact Orientation3 determinant sign for every finite binary64 input.
 *
 * Computes:
 *
 *     sign(det(b-a, c-a, d-a))
 *
 * Returns:
 *
 *     -1 negative
 *      0 coplanar
 *      1 positive
 *
 * Preconditions:
 *
 * - every input coordinate is finite.
 *
 * After binary64 coordinate decoding, all arithmetic is fixed-width integer
 * arithmetic.
 */
int orientationDyadicExact(
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


    const auto ax =
        decodeBinary64Coordinate(
            a.x
        );

    const auto ay =
        decodeBinary64Coordinate(
            a.y
        );

    const auto az =
        decodeBinary64Coordinate(
            a.z
        );


    const auto bx =
        decodeBinary64Coordinate(
            b.x
        );

    const auto by =
        decodeBinary64Coordinate(
            b.y
        );

    const auto bz =
        decodeBinary64Coordinate(
            b.z
        );


    const auto cx =
        decodeBinary64Coordinate(
            c.x
        );

    const auto cy =
        decodeBinary64Coordinate(
            c.y
        );

    const auto cz =
        decodeBinary64Coordinate(
            c.z
        );


    const auto dx =
        decodeBinary64Coordinate(
            d.x
        );

    const auto dy =
        decodeBinary64Coordinate(
            d.y
        );

    const auto dz =
        decodeBinary64Coordinate(
            d.z
        );


    const auto ux =
        subtractDyadicCoordinates(
            bx,
            ax
        );

    const auto uy =
        subtractDyadicCoordinates(
            by,
            ay
        );

    const auto uz =
        subtractDyadicCoordinates(
            bz,
            az
        );


    const auto vx =
        subtractDyadicCoordinates(
            cx,
            ax
        );

    const auto vy =
        subtractDyadicCoordinates(
            cy,
            ay
        );

    const auto vz =
        subtractDyadicCoordinates(
            cz,
            az
        );


    const auto wx =
        subtractDyadicCoordinates(
            dx,
            ax
        );

    const auto wy =
        subtractDyadicCoordinates(
            dy,
            ay
        );

    const auto wz =
        subtractDyadicCoordinates(
            dz,
            az
        );


    OrientationTripleMagnitude positive;
    OrientationTripleMagnitude negative;


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
    accumulateTriple(
        positive,
        negative,
        +1,
        ux, vy, wz
    );

    accumulateTriple(
        positive,
        negative,
        +1,
        uy, vz, wx
    );

    accumulateTriple(
        positive,
        negative,
        +1,
        uz, vx, wy
    );

    accumulateTriple(
        positive,
        negative,
        -1,
        uz, vy, wx
    );

    accumulateTriple(
        positive,
        negative,
        -1,
        uy, vx, wz
    );

    accumulateTriple(
        positive,
        negative,
        -1,
        ux, vz, wy
    );


    return compareUnsigned(
        positive,
        negative
    );
}


version(unittest)
{
    import std.bigint :
        BigInt;

    import std.bitmanip :
        DoubleRep;


    private alias P =
        Point3!double;


    /*
     * Independent exact finite binary64 decoder:
     *
     *     value = integer * 2^-1074
     */
    private BigInt oracleBinary64Integer(double value)
        @safe
    {
        assert(
            isFinite(value)
        );


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
     * Independent arbitrary-precision determinant oracle.
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


    private void checkCase(
        P a,
        P b,
        P c,
        P d
    )
        @safe
    {
        const int expected =
            oracleOrientation(
                a,
                b,
                c,
                d
            );

        const int actual =
            orientationDyadicExact(
                a,
                b,
                c,
                d
            );


        assert(
            actual ==
            expected
        );
    }


    private void checkExpected(
        int expected,
        P a,
        P b,
        P c,
        P d
    )
        @safe
    {
        assert(
            expected >= -1 &&
            expected <= 1
        );

        assert(
            oracleOrientation(
                a,
                b,
                c,
                d
            ) ==
            expected
        );

        assert(
            orientationDyadicExact(
                a,
                b,
                c,
                d
            ) ==
            expected
        );


        /*
         * One transposition is odd.
         */
        assert(
            orientationDyadicExact(
                a,
                c,
                b,
                d
            ) ==
            -expected
        );


        /*
         * A four-cycle is odd.
         */
        assert(
            orientationDyadicExact(
                b,
                c,
                d,
                a
            ) ==
            -expected
        );
    }


    private ulong nextRandom(ref ulong state)
        pure nothrow @safe @nogc
    {
        assert(
            state != 0
        );

        state ^=
            state << 13;

        state ^=
            state >> 7;

        state ^=
            state << 17;

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
                (
                    bits >>
                    52
                ) &
                0x7ffUL
            );


        /*
         * Raw exponent 2047 is NaN/infinity. Remap to largest finite
         * exponent while retaining fraction/sign entropy.
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


    private P randomFinitePoint(ref ulong state)
        pure nothrow @safe @nogc
    {
        return P(
            randomFiniteDouble(state),
            randomFiniteDouble(state),
            randomFiniteDouble(state)
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
     * Canonical positive geo3-d convention.
     */
    checkExpected(
        +1,
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


    /*
     * Exact coplanarity.
     */
    checkExpected(
        0,
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
            1.0,
            1.0,
            0.0
        )
    );


    /*
     * Complete degeneracy.
     */
    checkExpected(
        0,
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


    enum double minSubnormal =
        0x0.0000000000001p-1022;


    /*
     * Smallest positive subnormal along every basis direction.
     *
     * Ordinary binary64 triple multiplication cannot represent this
     * determinant; the dyadic backend must nevertheless return +1.
     */
    checkExpected(
        +1,
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
     * Determinant magnitude far beyond binary64.
     */
    checkExpected(
        +1,
        zero,
        P(
            double.max,
            0.0,
            0.0
        ),
        P(
            0.0,
            double.max,
            0.0
        ),
        P(
            0.0,
            0.0,
            double.max
        )
    );


    /*
     * Exact coordinate subtraction itself exceeds binary64.
     */
    checkExpected(
        +1,
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
     * Full-span exact coplanarity despite overflowing binary64 X/Y
     * differences.
     */
    checkExpected(
        0,
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
     * Extreme exponent imbalance.
     */
    checkExpected(
        +1,
        zero,
        P(
            double.max,
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
            double.max
        )
    );


    /*
     * Near-coplanar one-ULP positive case.
     */
    checkExpected(
        +1,
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
     * Corresponding representable value below the plane.
     */
    checkExpected(
        -1,
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
            0x1.fffffffffffffp+0
        )
    );
}


@safe unittest
{
    /*
     * The research prototype validates 5,000 arbitrary finite and 5,000
     * arbitrary full-range coplanar cases under DMD/LDC debug/release.
     *
     * Keep a deterministic production regression subset here.
     */
    alias P =
        Point3!double;

    ulong state =
        0xa54f_f53a_5f1d_36f1UL;


    enum randomCases =
        1_000;


    foreach (_; 0 .. randomCases)
    {
        checkCase(
            randomFinitePoint(state),
            randomFinitePoint(state),
            randomFinitePoint(state),
            randomFinitePoint(state)
        );
    }


    enum coplanarCases =
        1_000;


    foreach (_; 0 .. coplanarCases)
    {
        const P a =
            P(
                randomFiniteDouble(state),
                randomFiniteDouble(state),
                0.0
            );

        const P b =
            P(
                randomFiniteDouble(state),
                randomFiniteDouble(state),
                0.0
            );

        const P c =
            P(
                randomFiniteDouble(state),
                randomFiniteDouble(state),
                0.0
            );

        const P d =
            P(
                randomFiniteDouble(state),
                randomFiniteDouble(state),
                0.0
            );


        assert(
            oracleOrientation(
                a,
                b,
                c,
                d
            ) == 0
        );

        assert(
            orientationDyadicExact(
                a,
                b,
                c,
                d
            ) == 0
        );
    }
}
