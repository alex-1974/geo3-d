/**
 * Exact integral backend for three-dimensional orientation.
 *
 * INTERNAL IMPLEMENTATION MODULE.
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
module geo3.internal.orientation_integral;

import geo3.internal.fixed_uint :
    UIntFixed,
    addAssign,
    compareUnsigned,
    multiplyUnsigned;

import geo3.point :
    Point3;


/*
 * Exact signed coordinate difference represented as sign plus unsigned
 * fixed-width magnitude.
 */
private struct SignedDiff(size_t Limbs)
if (Limbs > 0)
{
    int sign;
    UIntFixed!Limbs magnitude;
}


/*
 * Exact unsigned magnitude without overflowing int.min.
 */
private uint unsignedMagnitude(int value)
    pure nothrow @safe @nogc
{
    if (value >= 0)
        return cast(uint) value;

    return
        cast(uint)(
            -(value + 1)
        ) +
        1U;
}


/*
 * Exact unsigned magnitude without overflowing long.min.
 */
private ulong unsignedMagnitude(long value)
    pure nothrow @safe @nogc
{
    if (value >= 0)
        return cast(ulong) value;

    return
        cast(ulong)(
            -(value + 1)
        ) +
        1UL;
}


/*
 * Exact lhs-rhs for int coordinates.
 *
 * The complete difference magnitude is at most 2^32-1 and therefore fits in
 * one uint limb.
 */
private SignedDiff!1 signedDifference(
    int lhs,
    int rhs
)
    pure nothrow @safe @nogc
{
    SignedDiff!1 result;

    if (lhs == rhs)
        return result;

    result.sign =
        lhs > rhs
            ? 1
            : -1;

    uint magnitude;

    if ((lhs < 0) != (rhs < 0))
    {
        magnitude =
            unsignedMagnitude(lhs) +
            unsignedMagnitude(rhs);
    }
    else if (lhs > rhs)
    {
        magnitude =
            cast(uint)(
                lhs - rhs
            );
    }
    else
    {
        magnitude =
            cast(uint)(
                rhs - lhs
            );
    }

    result.magnitude.limb[0] =
        magnitude;

    return result;
}


/*
 * Exact lhs-rhs for long coordinates.
 *
 * The complete difference magnitude is at most 2^64-1 and therefore fits in
 * two uint limbs.
 */
private SignedDiff!2 signedDifference(
    long lhs,
    long rhs
)
    pure nothrow @safe @nogc
{
    SignedDiff!2 result;

    if (lhs == rhs)
        return result;

    result.sign =
        lhs > rhs
            ? 1
            : -1;

    ulong magnitude;

    if ((lhs < 0) != (rhs < 0))
    {
        magnitude =
            unsignedMagnitude(lhs) +
            unsignedMagnitude(rhs);
    }
    else if (lhs > rhs)
    {
        magnitude =
            cast(ulong)(
                lhs - rhs
            );
    }
    else
    {
        magnitude =
            cast(ulong)(
                rhs - lhs
            );
    }

    result.magnitude.limb[0] =
        cast(uint) magnitude;

    result.magnitude.limb[1] =
        cast(uint)(
            magnitude >> 32
        );

    return result;
}


/*
 * Adds one exact signed determinant term
 *
 *     coefficientSign * x * y * z
 *
 * to the corresponding unsigned sign bucket.
 */
private void accumulateTerm(
    size_t BucketLimbs,
    size_t DiffLimbs
)(
    ref UIntFixed!BucketLimbs positive,
    ref UIntFixed!BucketLimbs negative,
    int coefficientSign,
    ref const SignedDiff!DiffLimbs x,
    ref const SignedDiff!DiffLimbs y,
    ref const SignedDiff!DiffLimbs z
)
    pure nothrow @safe @nogc
if (BucketLimbs >= 3 * DiffLimbs)
{
    assert(
        coefficientSign == -1 ||
        coefficientSign == 1
    );

    if (
        x.sign == 0 ||
        y.sign == 0 ||
        z.sign == 0
    )
        return;

    const int sign =
        coefficientSign *
        x.sign *
        y.sign *
        z.sign;

    const auto xy =
        multiplyUnsigned(
            x.magnitude,
            y.magnitude
        );

    const auto xyz =
        multiplyUnsigned(
            xy,
            z.magnitude
        );

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
     * Width proofs:
     *
     * int:
     *   each triple product < 2^96;
     *   at most five non-zero terms share one sign;
     *   bucket < 5 * 2^96 < 2^99;
     *   four uint limbs provide 128 bits.
     *
     * long:
     *   each triple product < 2^192;
     *   bucket < 5 * 2^192 < 2^195;
     *   seven uint limbs provide 224 bits.
     */
    assert(success);
}


/*
 * Exact determinant sign for:
 *
 *     det(b-a, c-a, d-a)
 *
 * using positive and negative unsigned Leibniz-term buckets.
 */
private int orientationIntegralSignImpl(
    T,
    size_t DiffLimbs,
    size_t BucketLimbs
)(
    Point3!T a,
    Point3!T b,
    Point3!T c,
    Point3!T d
)
    pure nothrow @safe @nogc
if (
    (is(T == int) || is(T == long)) &&
    BucketLimbs >= 3 * DiffLimbs
)
{
    const auto ux =
        signedDifference(
            b.x,
            a.x
        );

    const auto uy =
        signedDifference(
            b.y,
            a.y
        );

    const auto uz =
        signedDifference(
            b.z,
            a.z
        );


    const auto vx =
        signedDifference(
            c.x,
            a.x
        );

    const auto vy =
        signedDifference(
            c.y,
            a.y
        );

    const auto vz =
        signedDifference(
            c.z,
            a.z
        );


    const auto wx =
        signedDifference(
            d.x,
            a.x
        );

    const auto wy =
        signedDifference(
            d.y,
            a.y
        );

    const auto wz =
        signedDifference(
            d.z,
            a.z
        );


    UIntFixed!BucketLimbs positive;
    UIntFixed!BucketLimbs negative;


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
    accumulateTerm(
        positive,
        negative,
        +1,
        ux, vy, wz
    );

    accumulateTerm(
        positive,
        negative,
        +1,
        uy, vz, wx
    );

    accumulateTerm(
        positive,
        negative,
        +1,
        uz, vx, wy
    );

    accumulateTerm(
        positive,
        negative,
        -1,
        uz, vy, wx
    );

    accumulateTerm(
        positive,
        negative,
        -1,
        uy, vx, wz
    );

    accumulateTerm(
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


/**
 * Exact orient3d determinant sign over the complete int coordinate domain.
 */
int orientationIntegralSign(
    Point3!int a,
    Point3!int b,
    Point3!int c,
    Point3!int d
)
    pure nothrow @safe @nogc
{
    return orientationIntegralSignImpl!(
        int,
        1,
        4
    )(
        a,
        b,
        c,
        d
    );
}


/**
 * Exact orient3d determinant sign over the complete long coordinate domain.
 */
int orientationIntegralSign(
    Point3!long a,
    Point3!long b,
    Point3!long c,
    Point3!long d
)
    pure nothrow @safe @nogc
{
    return orientationIntegralSignImpl!(
        long,
        2,
        7
    )(
        a,
        b,
        c,
        d
    );
}
