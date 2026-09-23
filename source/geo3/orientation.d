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
 * Floating-point orientation is deliberately not introduced by the integral
 * implementation slice.
 *
 * It will become public only after the certified filter and exact fallback
 * pipeline is present.
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
    !__traits(
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
