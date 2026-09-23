module geo.internal.fixed_uint;


/*
 * INTERNAL IMPLEMENTATION MODULE.
 *
 * Small fixed-width unsigned integer primitives for robust geometric
 * predicates and constructions.
 *
 * This is deliberately not a general arbitrary-precision integer type.
 *
 * Storage is:
 *
 * - fixed at compile time;
 * - inline/value-owned;
 * - little-endian in base 2^32;
 * - allocation-free.
 */


/**
 * Fixed-width unsigned integer using little-endian 32-bit limbs.
 */
struct UIntFixed(size_t Limbs)
if (Limbs > 0)
{
    uint[Limbs] limb;


    @property bool isZero() const
        pure nothrow @safe @nogc
    {
        foreach (value; limb)
        {
            if (value != 0)
                return false;
        }

        return true;
    }
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

    while (index > 0)
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
 * Exact fixed-width unsigned addition.
 *
 * The caller guarantees that the mathematical result fits in the
 * selected width.
 */
UIntFixed!Limbs addUnsigned(size_t Limbs)(
    ref const UIntFixed!Limbs lhs,
    ref const UIntFixed!Limbs rhs
)
    pure nothrow @safe @nogc
{
    /*
     * Begin with lhs verbatim and touch only the active span of rhs.
     *
     * This avoids traversing large zero ranges in fixed-width values such
     * as exact dyadic products.
     */
    UIntFixed!Limbs result =
        lhs;

    const size_t rhsFirst =
        firstNonZeroLimb(rhs);

    if (rhsFirst == Limbs)
        return result;

    const size_t rhsEnd =
        pastLastNonZeroLimb(rhs);

    ulong carry = 0;

    foreach (index; rhsFirst .. rhsEnd)
    {
        const ulong sum =
            cast(ulong) result.limb[index] +
            cast(ulong) rhs.limb[index] +
            carry;

        result.limb[index] =
            cast(uint) sum;

        carry =
            sum >> 32;
    }

    size_t index =
        rhsEnd;

    while (carry != 0)
    {
        assert(index < Limbs);

        const ulong sum =
            cast(ulong) result.limb[index] +
            carry;

        result.limb[index] =
            cast(uint) sum;

        carry =
            sum >> 32;

        ++index;
    }

    return result;
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
        compareUnsigned(lhs, rhs) >= 0
    );

    /*
     * Preserve every limb outside rhs's active range directly from lhs.
     *
     * Below the least-significant non-zero rhs limb no borrow can exist.
     * Above the active range only a remaining borrow needs propagation.
     */
    UIntFixed!Limbs result =
        lhs;

    const size_t rhsFirst =
        firstNonZeroLimb(rhs);

    if (rhsFirst == Limbs)
        return result;

    const size_t rhsEnd =
        pastLastNonZeroLimb(rhs);

    ulong borrow = 0;

    foreach (index; rhsFirst .. rhsEnd)
    {
        const ulong lhsValue =
            cast(ulong) result.limb[index];

        const ulong rhsValue =
            cast(ulong) rhs.limb[index] +
            borrow;

        if (lhsValue >= rhsValue)
        {
            result.limb[index] =
                cast(uint)(
                    lhsValue - rhsValue
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

    size_t index =
        rhsEnd;

    while (borrow != 0)
    {
        assert(index < Limbs);

        if (result.limb[index] != 0)
        {
            --result.limb[index];
            borrow = 0;
        }
        else
        {
            result.limb[index] =
                uint.max;

            ++index;
        }
    }

    return result;
}


/*
 * Index of the least-significant non-zero limb.
 *
 * Returns Limbs for zero.
 */
private size_t firstNonZeroLimb(size_t Limbs)(
    ref const UIntFixed!Limbs value
)
    pure nothrow @safe @nogc
{
    foreach (i; 0 .. Limbs)
    {
        if (value.limb[i] != 0)
            return i;
    }

    return Limbs;
}


/*
 * One-past the most-significant non-zero limb.
 *
 * Returns zero for zero.
 */
private size_t pastLastNonZeroLimb(size_t Limbs)(
    ref const UIntFixed!Limbs value
)
    pure nothrow @safe @nogc
{
    for (size_t i = Limbs; i != 0; --i)
    {
        if (value.limb[i - 1] != 0)
            return i;
    }

    return 0;
}


/**
 * Exact fixed-width multiplication.
 *
 * Multiplying an M-limb value by an N-limb value yields an
 * (M + N)-limb result.
 *
 * The implementation uses ordinary base-2^32 schoolbook
 * multiplication.
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


    /*
     * Fixed-width storage is intentionally sized for the complete
     * numerical domain, but ordinary dyadic values usually occupy only
     * a small contiguous range of limbs.
     *
     * Skip leading and trailing zero ranges while retaining every limb
     * inside the active span. Interior zero limbs must still be
     * processed because carry may propagate through them.
     */
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
             * Maximum possible accumulator:
             *
             *     (2^32 - 1)^2
             *   + (2^32 - 1)
             *   + (2^32 - 1)
             *
             * = 2^64 - 1
             *
             * therefore ulong is exactly sufficient.
             */
            const ulong accumulated =
                cast(ulong) lhsWord *
                    cast(ulong) rhs.limb[j]
                + cast(ulong)
                    result.limb[index]
                + carry;

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
         * No earlier multiplication row can have reached this limb.
         */
        assert(
            result.limb[carryIndex] == 0
        );

        result.limb[carryIndex] =
            cast(uint) carry;
    }

    return result;
}


@safe unittest
{
    /*
     * Sparse operands near the high end of their fixed-width storage.
     */
    {
        UIntFixed!4 a;
        UIntFixed!5 b;

        a.limb[3] = 2;
        b.limb[4] = 3;

        const auto product =
            multiplyUnsigned(a, b);

        foreach (i; 0 .. 7)
            assert(product.limb[i] == 0);

        assert(product.limb[7] == 6);
        assert(product.limb[8] == 0);
    }


    /*
     * Carry must propagate through an interior zero limb.
     *
     * Skipping zero limbs inside the active range would make this
     * result incorrect.
     */
    {
        UIntFixed!1 a;
        UIntFixed!3 b;

        a.limb[0] = uint.max;

        b.limb[0] = uint.max;
        b.limb[1] = 0;
        b.limb[2] = 1;

        const auto product =
            multiplyUnsigned(a, b);

        assert(product.limb[0] == 1);
        assert(
            product.limb[1] ==
            uint.max - 1
        );
        assert(
            product.limb[2] ==
            uint.max
        );
        assert(product.limb[3] == 0);
    }


    /*
     * Sparse multiplication preserves unequal-width positioning.
     */
    {
        UIntFixed!6 a;
        UIntFixed!2 b;

        a.limb[4] = 7;
        b.limb[1] = 9;

        const auto product =
            multiplyUnsigned(a, b);

        assert(product.limb[5] == 63);

        foreach (i; 0 .. 5)
            assert(product.limb[i] == 0);

        foreach (i; 6 .. 8)
            assert(product.limb[i] == 0);
    }


    /*
     * Initial value and zero detection.
     */
    {
        UIntFixed!2 value;

        assert(value.isZero);

        value.limb[0] = 1;

        assert(!value.isZero);
    }


    /*
     * Comparison proceeds from the most significant limb.
     */
    {
        UIntFixed!2 a;
        UIntFixed!2 b;

        a.limb[0] = uint.max;

        b.limb[1] = 1;

        assert(
            compareUnsigned(a, b) < 0
        );

        assert(
            compareUnsigned(b, a) > 0
        );

        assert(
            compareUnsigned(a, a) == 0
        );
    }


    /*
     * Addition propagates carry across limbs.
     */
    {
        UIntFixed!2 a;
        UIntFixed!2 b;

        a.limb[0] =
            uint.max;

        b.limb[0] = 1;

        const auto sum =
            addUnsigned(a, b);

        assert(sum.limb[0] == 0);
        assert(sum.limb[1] == 1);
    }


    /*
     * Subtraction propagates borrow across limbs.
     */
    {
        UIntFixed!2 a;
        UIntFixed!2 b;

        a.limb[1] = 1;
        b.limb[0] = 1;

        const auto difference =
            subtractUnsigned(a, b);

        assert(
            difference.limb[0] ==
            uint.max
        );

        assert(
            difference.limb[1] == 0
        );
    }


    /*
     * 1-limb multiplication.
     */
    {
        UIntFixed!1 a;
        UIntFixed!1 b;

        a.limb[0] =
            uint.max;

        b.limb[0] =
            uint.max;

        const auto product =
            multiplyUnsigned(a, b);

        static assert(
            is(
                typeof(
                    multiplyUnsigned(a, b)
                ) ==
                UIntFixed!2
            )
        );

        assert(
            product.limb[0] == 1
        );

        assert(
            product.limb[1] ==
            uint.max - 1
        );
    }


    /*
     * Multiplication works across unequal widths.
     *
     *     2^32 * 1
     *   = 2^32
     */
    {
        UIntFixed!2 a;
        UIntFixed!1 b;

        a.limb[1] = 1;
        b.limb[0] = 1;

        const auto product =
            multiplyUnsigned(a, b);

        static assert(
            is(
                typeof(
                    multiplyUnsigned(a, b)
                ) ==
                UIntFixed!3
            )
        );

        assert(product.limb[0] == 0);
        assert(product.limb[1] == 1);
        assert(product.limb[2] == 0);
    }


    /*
     * Multi-limb carry propagation.
     *
     *     (2^64 - 1) * (2^32 - 1)
     */
    {
        UIntFixed!2 a;
        UIntFixed!1 b;

        a.limb[0] =
            uint.max;

        a.limb[1] =
            uint.max;

        b.limb[0] =
            uint.max;

        const auto product =
            multiplyUnsigned(a, b);

        assert(
            product.limb[0] == 1
        );

        assert(
            product.limb[1] ==
            uint.max
        );

        assert(
            product.limb[2] ==
            uint.max - 1
        );
    }
}


@safe unittest
{
    /*
     * Sparse addition must propagate carry beyond rhs's active span.
     */
    {
        UIntFixed!6 lhs;
        UIntFixed!6 rhs;

        lhs.limb[2] = uint.max;
        lhs.limb[3] = uint.max;

        rhs.limb[2] = 1;

        const auto result =
            addUnsigned(
                lhs,
                rhs
            );

        assert(result.limb[0] == 0);
        assert(result.limb[1] == 0);
        assert(result.limb[2] == 0);
        assert(result.limb[3] == 0);
        assert(result.limb[4] == 1);
        assert(result.limb[5] == 0);
    }


    /*
     * Sparse subtraction must propagate borrow beyond rhs's active span.
     */
    {
        UIntFixed!6 lhs;
        UIntFixed!6 rhs;

        lhs.limb[4] = 1;
        rhs.limb[2] = 1;

        const auto result =
            subtractUnsigned(
                lhs,
                rhs
            );

        assert(result.limb[0] == 0);
        assert(result.limb[1] == 0);
        assert(result.limb[2] == uint.max);
        assert(result.limb[3] == uint.max);
        assert(result.limb[4] == 0);
        assert(result.limb[5] == 0);
    }
}
