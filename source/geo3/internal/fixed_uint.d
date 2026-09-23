/**
 * Small fixed-width unsigned integer primitives for robust 3D predicates.
 *
 * This is deliberately not a general arbitrary-precision integer type.
 *
 * Storage is:
 *
 * - fixed at compile time;
 * - inline/value-owned;
 * - little-endian in base 2^32;
 * - allocation-free.
 *
 * Portions of the larger-width arithmetic are derived from geo-d v2.0.0:
 *
 *     repository: https://github.com/alex-1974/geo-d
 *     tag:        v2.0.0
 *     commit:     83c974e0ca018bb378ce65a08ef9d3178decf9ad
 *     source:     source/geo/internal/fixed_uint.d
 *     blob:       d50a7dc640584e3a4a95d6c49efa2007616b8104
 *
 * No production dependency on geo-d is introduced.
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
module geo3.internal.fixed_uint;


/**
 * Fixed-width unsigned integer using little-endian 32-bit limbs.
 */
struct UIntFixed(size_t Limbs)
if (Limbs > 0)
{
    uint[Limbs] limb;
}


/**
 * Exact unsigned comparison.
 *
 * Returns:
 *
 *     -1 lhs < rhs
 *      0 lhs == rhs
 *      1 lhs > rhs
 */
int compareUnsigned(size_t Limbs)(
    ref const UIntFixed!Limbs lhs,
    ref const UIntFixed!Limbs rhs
)
    pure nothrow @safe @nogc
{
    size_t index =
        Limbs;

    while (index != 0)
    {
        --index;

        if (lhs.limb[index] < rhs.limb[index])
            return -1;

        if (lhs.limb[index] > rhs.limb[index])
            return 1;
    }

    return 0;
}


/**
 * Adds an unsigned fixed-width value into a wider or equally wide target.
 *
 * Returns false only when the mathematical result does not fit in the
 * target width.
 */
bool addAssign(
    size_t TargetLimbs,
    size_t ValueLimbs
)(
    ref UIntFixed!TargetLimbs target,
    ref const UIntFixed!ValueLimbs value
)
    pure nothrow @safe @nogc
if (TargetLimbs >= ValueLimbs)
{
    ulong carry = 0;

    foreach (index; 0 .. ValueLimbs)
    {
        const ulong sum =
            cast(ulong) target.limb[index] +
            cast(ulong) value.limb[index] +
            carry;

        target.limb[index] =
            cast(uint) sum;

        carry =
            sum >> 32;
    }


    size_t index =
        ValueLimbs;

    while (
        carry != 0 &&
        index < TargetLimbs
    )
    {
        const ulong sum =
            cast(ulong) target.limb[index] +
            carry;

        target.limb[index] =
            cast(uint) sum;

        carry =
            sum >> 32;

        ++index;
    }


    return carry == 0;
}


/**
 * Exact fixed-width unsigned subtraction.
 *
 * Precondition:
 *
 *     lhs >= rhs
 */
UIntFixed!Limbs subtractUnsigned(size_t Limbs)(
    ref const UIntFixed!Limbs lhs,
    ref const UIntFixed!Limbs rhs
)
    pure nothrow @safe @nogc
{
    assert(
        compareUnsigned(
            lhs,
            rhs
        ) >= 0
    );


    UIntFixed!Limbs result;
    ulong borrow = 0;


    foreach (index; 0 .. Limbs)
    {
        const ulong lhsValue =
            cast(ulong) lhs.limb[index];

        const ulong rhsValue =
            cast(ulong) rhs.limb[index] +
            borrow;


        if (lhsValue >= rhsValue)
        {
            result.limb[index] =
                cast(uint)(
                    lhsValue -
                    rhsValue
                );

            borrow = 0;
        }
        else
        {
            result.limb[index] =
                cast(uint)(
                    0x1_0000_0000UL +
                    lhsValue -
                    rhsValue
                );

            borrow = 1;
        }
    }


    assert(borrow == 0);

    return result;
}


/*
 * Index of the least-significant non-zero limb.
 *
 * Returns Limbs for exact zero.
 */
private size_t firstNonZeroLimb(size_t Limbs)(
    ref const UIntFixed!Limbs value
)
    pure nothrow @safe @nogc
{
    foreach (index; 0 .. Limbs)
    {
        if (value.limb[index] != 0)
            return index;
    }

    return Limbs;
}


/*
 * One-past the most-significant non-zero limb.
 *
 * Returns zero for exact zero.
 */
private size_t pastLastNonZeroLimb(size_t Limbs)(
    ref const UIntFixed!Limbs value
)
    pure nothrow @safe @nogc
{
    size_t index =
        Limbs;

    while (index != 0)
    {
        if (value.limb[index - 1] != 0)
            return index;

        --index;
    }

    return 0;
}


/**
 * Exact unsigned fixed-width multiplication.
 *
 * An M-limb value multiplied by an N-limb value produces an exact
 * (M + N)-limb result.
 *
 * Leading and trailing zero ranges are skipped. Interior zero limbs remain
 * part of the arithmetic because carry can propagate through them.
 */
UIntFixed!(LhsLimbs + RhsLimbs) multiplyUnsigned(
    size_t LhsLimbs,
    size_t RhsLimbs
)(
    ref const UIntFixed!LhsLimbs lhs,
    ref const UIntFixed!RhsLimbs rhs
)
    pure nothrow @safe @nogc
{
    UIntFixed!(LhsLimbs + RhsLimbs) result;


    const size_t lhsFirst =
        firstNonZeroLimb(lhs);

    if (lhsFirst == LhsLimbs)
        return result;


    const size_t rhsFirst =
        firstNonZeroLimb(rhs);

    if (rhsFirst == RhsLimbs)
        return result;


    const size_t lhsEnd =
        pastLastNonZeroLimb(lhs);

    const size_t rhsEnd =
        pastLastNonZeroLimb(rhs);


    foreach (i; lhsFirst .. lhsEnd)
    {
        const uint lhsWord =
            lhs.limb[i];

        if (lhsWord == 0)
            continue;


        ulong carry = 0;


        foreach (j; rhsFirst .. rhsEnd)
        {
            const size_t index =
                i + j;

            /*
             * Maximum:
             *
             *     (2^32 - 1)^2
             *       + (2^32 - 1)
             *       + (2^32 - 1)
             *     = 2^64 - 1
             */
            const ulong accumulated =
                cast(ulong) lhsWord *
                    cast(ulong) rhs.limb[j] +
                cast(ulong) result.limb[index] +
                carry;

            result.limb[index] =
                cast(uint) accumulated;

            carry =
                accumulated >> 32;
        }


        const size_t carryIndex =
            i + rhsEnd;

        assert(
            carryIndex <
            LhsLimbs + RhsLimbs
        );

        /*
         * An earlier multiplication row cannot yet have reached this limb.
         */
        assert(
            result.limb[carryIndex] ==
            0
        );

        result.limb[carryIndex] =
            cast(uint) carry;
    }


    return result;
}


@safe unittest
{
    UIntFixed!1 one;
    one.limb[0] = 1;

    UIntFixed!1 max32;
    max32.limb[0] = uint.max;


    const auto product =
        multiplyUnsigned(
            max32,
            max32
        );

    assert(
        product.limb[0] ==
        1
    );

    assert(
        product.limb[1] ==
        uint.max - 1
    );


    UIntFixed!2 sum;

    assert(
        addAssign(
            sum,
            max32
        )
    );

    assert(
        addAssign(
            sum,
            one
        )
    );

    assert(
        sum.limb[0] ==
        0
    );

    assert(
        sum.limb[1] ==
        1
    );


    UIntFixed!2 smaller;
    smaller.limb[0] =
        uint.max;

    assert(
        compareUnsigned(
            sum,
            smaller
        ) > 0
    );


    /*
     * Borrow across a complete zero limb.
     *
     *     2^64 - 1
     */
    UIntFixed!3 large;
    UIntFixed!3 subtract;

    large.limb[2] =
        1;

    subtract.limb[0] =
        1;

    const auto difference =
        subtractUnsigned(
            large,
            subtract
        );

    assert(
        difference.limb[0] ==
        uint.max
    );

    assert(
        difference.limb[1] ==
        uint.max
    );

    assert(
        difference.limb[2] ==
        0
    );


    /*
     * Sparse high-limb multiplication.
     */
    UIntFixed!4 sparseA;
    UIntFixed!5 sparseB;

    sparseA.limb[3] =
        2;

    sparseB.limb[4] =
        3;

    const auto sparseProduct =
        multiplyUnsigned(
            sparseA,
            sparseB
        );

    assert(
        sparseProduct.limb[7] ==
        6
    );

    foreach (index; 0 .. 7)
    {
        assert(
            sparseProduct.limb[index] ==
            0
        );
    }

    assert(
        sparseProduct.limb[8] ==
        0
    );


    /*
     * Carry through an interior zero limb.
     */
    UIntFixed!1 carryA;
    UIntFixed!3 carryB;

    carryA.limb[0] =
        uint.max;

    carryB.limb[0] =
        uint.max;

    carryB.limb[1] =
        0;

    carryB.limb[2] =
        1;

    const auto carryProduct =
        multiplyUnsigned(
            carryA,
            carryB
        );

    assert(
        carryProduct.limb[0] ==
        1
    );

    assert(
        carryProduct.limb[1] ==
        uint.max - 1
    );

    assert(
        carryProduct.limb[2] ==
        uint.max
    );

    assert(
        carryProduct.limb[3] ==
        0
    );
}
