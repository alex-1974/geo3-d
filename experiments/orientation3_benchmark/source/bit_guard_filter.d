module bit_guard_filter;

import geo3.internal.binary64_rounding :
    roundedAdd,
    roundedMul,
    roundedSub;

import geo3.internal.orientation_filter :
    OrientationFilterResult;

import geo3.point :
    Point3;

import std.bitmanip :
    DoubleRep;


/*
 * EXPERIMENT ONLY.
 *
 * Bit-level binary64 classification candidates for the production
 * Orientation3 filter.
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
        /*
         * exponent == 0:
         *
         * fraction == 0 -> +/- zero
         * fraction != 0 -> subnormal
         */
        return
            representation.fraction ==
            0;
    }


    /*
     * 1 .. 2046 -> normal
     * 2047      -> infinity / NaN
     */
    return
        exponent !=
        0x7ff;
}


/*
 * Valid multiplication result for the certified filter.
 *
 * Normal products are accepted.
 *
 * Zero is accepted only if at least one operand is mathematically zero.
 * Thus a product that rounded to zero from two non-zero operands is rejected
 * as underflow.
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
 * The error bound must be a finite normal non-zero binary64 value.
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


pragma(inline, true)
private double absolute(double value)
    pure nothrow @safe @nogc
{
    return value < 0.0
        ? -value
        : value;
}


/*
 * Candidate preserving the production filter's non-finite-input behaviour.
 */
OrientationFilterResult bitGuardFilter(
    Point3!double a,
    Point3!double b,
    Point3!double c,
    Point3!double d
)
    pure nothrow @safe @nogc
{
    if (
        !finiteBits(a.x) ||
        !finiteBits(a.y) ||
        !finiteBits(a.z) ||
        !finiteBits(b.x) ||
        !finiteBits(b.y) ||
        !finiteBits(b.z) ||
        !finiteBits(c.x) ||
        !finiteBits(c.y) ||
        !finiteBits(c.z) ||
        !finiteBits(d.x) ||
        !finiteBits(d.y) ||
        !finiteBits(d.z)
    )
    {
        return
            OrientationFilterResult.uncertain;
    }


    return bitGuardFilterFinite(
        a,
        b,
        c,
        d
    );
}


/*
 * Candidate for callers that already satisfy the finite-coordinate
 * precondition.
 *
 * Production orientation(Point3!double, ...) already declares that
 * precondition. This variant lets the benchmark quantify the cost of
 * rechecking it inside the filter.
 */
OrientationFilterResult bitGuardFilterFinite(
    Point3!double a,
    Point3!double b,
    Point3!double c,
    Point3!double d
)
    pure nothrow @safe @nogc
{
    enum double errboundA =
        0x1.c000000000007p-51;


    const double adx =
        roundedSub(a.x, d.x);

    const double bdx =
        roundedSub(b.x, d.x);

    const double cdx =
        roundedSub(c.x, d.x);

    const double ady =
        roundedSub(a.y, d.y);

    const double bdy =
        roundedSub(b.y, d.y);

    const double cdy =
        roundedSub(c.y, d.y);

    const double adz =
        roundedSub(a.z, d.z);

    const double bdz =
        roundedSub(b.z, d.z);

    const double cdz =
        roundedSub(c.z, d.z);


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
        roundedMul(bdx, cdy);

    const double cdxbdy =
        roundedMul(cdx, bdy);

    const double cdxady =
        roundedMul(cdx, ady);

    const double adxcdy =
        roundedMul(adx, cdy);

    const double adxbdy =
        roundedMul(adx, bdy);

    const double bdxady =
        roundedMul(bdx, ady);


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
            errboundA,
            permanent
        );


    if (!validErrorBoundBits(errbound))
    {
        return
            OrientationFilterResult.uncertain;
    }


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
