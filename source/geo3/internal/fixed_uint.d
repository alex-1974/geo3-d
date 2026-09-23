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
    size_t index = Limbs;

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
 * Exact unsigned fixed-width multiplication.
 *
 * The result has exactly the sum of the operand limb counts and therefore
 * always has enough storage for the full mathematical product.
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

    foreach (i; 0 .. LhsLimbs)
    {
        ulong carry = 0;

        foreach (j; 0 .. RhsLimbs)
        {
            const size_t index =
                i + j;

            /*
             * Maximum possible value is exactly <= ulong.max:
             *
             *     (2^32 - 1)^2
             *       + (2^32 - 1)
             *       + (2^32 - 1)
             *     = 2^64 - 1
             */
            const ulong accumulated =
                cast(ulong) lhs.limb[i] *
                    cast(ulong) rhs.limb[j] +
                cast(ulong) result.limb[index] +
                carry;

            result.limb[index] =
                cast(uint) accumulated;

            carry =
                accumulated >> 32;
        }

        size_t index =
            i + RhsLimbs;

        while (carry != 0)
        {
            assert(
                index <
                LhsLimbs + RhsLimbs
            );

            const ulong accumulated =
                cast(ulong) result.limb[index] +
                carry;

            result.limb[index] =
                cast(uint) accumulated;

            carry =
                accumulated >> 32;

            ++index;
        }
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

    /*
     * (2^32 - 1)^2
     *
     *     = (2^32 - 2) * 2^32 + 1
     */
    assert(product.limb[0] == 1);
    assert(product.limb[1] == uint.max - 1);

    UIntFixed!2 sum;

    assert(addAssign(sum, max32));
    assert(addAssign(sum, one));

    assert(sum.limb[0] == 0);
    assert(sum.limb[1] == 1);

    UIntFixed!2 smaller;
    smaller.limb[0] = uint.max;

    assert(
        compareUnsigned(sum, smaller) >
        0
    );
}
