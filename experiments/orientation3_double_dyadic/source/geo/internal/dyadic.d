module geo.internal.dyadic;

import geo.internal.fixed_uint :
    UIntFixed,
    addUnsigned,
    compareUnsigned,
    multiplyUnsigned,
    subtractUnsigned;

import std.bitmanip :
    DoubleRep;

import std.math.traits :
    isFinite;


/*
 * INTERNAL IMPLEMENTATION MODULE.
 *
 * Exact dyadic representation of finite binary64 coordinates.
 *
 * Every finite IEEE binary64 value can be represented exactly as:
 *
 *     integer * 2^-1074
 *
 * The integer magnitude therefore needs at most 2098 bits.
 *
 * 66 x 32 = 2112 bits, which is sufficient for every finite double.
 *
 * This module provides representation and decoding only. It deliberately
 * does not implement general signed arithmetic or rational arithmetic.
 */


enum size_t dyadicCoordinateLimbs = 66;


/*
 * Every exact coordinate is represented in units of 2^-1074.
 *
 * Therefore an integral coordinate n is embedded as:
 *
 *     n << 1074
 */
enum uint dyadicScaleShift = 1074;


alias DyadicCoordinateMagnitude =
    UIntFixed!dyadicCoordinateLimbs;


/**
 * Exact signed dyadic coordinate represented in units of 2^-1074.
 *
 * sign:
 *
 *     -1 negative
 *      0 zero
 *      1 positive
 */
struct SignedDyadicCoordinate
{
    int sign;
    DyadicCoordinateMagnitude magnitude;
}


/*
 * Places one non-zero unsigned 64-bit value at the requested bit
 * position.
 *
 * The destination is initially zero, so bitwise OR is sufficient.
 *
 * The caller is responsible for choosing a shift for which the value
 * fits in DyadicCoordinateMagnitude. Boundary assertions below guard
 * the final limbs.
 */
private void setShiftedUnsigned64(
    ref DyadicCoordinateMagnitude result,
    ulong mantissa,
    uint shift
)
    pure nothrow @safe @nogc
{
    assert(mantissa != 0);
    assert(shift <= 2045);

    const size_t base =
        shift / 32;

    const uint offset =
        shift % 32;

    const uint lower =
        cast(uint) mantissa;

    const uint upper =
        cast(uint)(mantissa >> 32);

    const ulong shiftedLower =
        cast(ulong) lower << offset;

    result.limb[base] |=
        cast(uint) shiftedLower;

    if (base + 1 < dyadicCoordinateLimbs)
    {
        result.limb[base + 1] |=
            cast(uint)(
                shiftedLower >> 32
            );
    }

    if (upper != 0)
    {
        const ulong shiftedUpper =
            cast(ulong) upper << offset;

        assert(
            base + 1 <
            dyadicCoordinateLimbs
        );

        result.limb[base + 1] |=
            cast(uint) shiftedUpper;

        if (
            base + 2 <
            dyadicCoordinateLimbs
        )
        {
            result.limb[base + 2] |=
                cast(uint)(
                    shiftedUpper >> 32
                );
        }
        else
        {
            assert(
                (shiftedUpper >> 32) == 0
            );
        }
    }
}


/**
 * Decodes one finite binary64 value exactly as a signed integer
 * multiple of 2^-1074.
 *
 * For subnormals:
 *
 *     value = fraction * 2^-1074
 *
 * For normals:
 *
 *     value =
 *         (2^52 + fraction)
 *         * 2^(rawExponent - 1075)
 *
 * Therefore, in units of 2^-1074:
 *
 *     integer =
 *         mantissa << (rawExponent - 1)
 *
 * Positive and negative zero both produce sign == 0.
 */
SignedDyadicCoordinate decodeBinary64Coordinate(
    double value
)
    pure nothrow @safe @nogc
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
        /*
         * Zero or subnormal.
         */
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

    SignedDyadicCoordinate result;

    if (mantissa == 0)
    {
        /*
         * +0 and -0 are the same mathematical coordinate.
         */
        result.sign = 0;
        return result;
    }

    result.sign =
        representation.sign
            ? -1
            : 1;

    setShiftedUnsigned64(
        result.magnitude,
        mantissa,
        shift
    );

    return result;
}


/*
 * Exact unsigned magnitude without overflowing the signed source type.
 */
private ulong integralMagnitude(T)(
    T value
)
    pure nothrow @safe @nogc
if (is(T == int) || is(T == long))
{
    if (value >= 0)
        return cast(ulong) value;

    static if (is(T == int))
    {
        /*
         * int fits completely in long, so negation after widening is
         * safe even for int.min.
         */
        return cast(ulong)(
            -cast(long) value
        );
    }
    else
    {
        /*
         * Negating long.min directly would overflow.
         *
         * For every negative long:
         *
         *     |value| = -(value + 1) + 1
         *
         * The signed negation is now representable and the final +1 is
         * performed as ulong.
         */
        return
            cast(ulong)(
                -(value + 1)
            ) +
            1UL;
    }
}


/*
 * Exact integral coordinate embedding into the common 2^-1074 scale.
 */
private SignedDyadicCoordinate decodeIntegralCoordinate(T)(
    T value
)
    pure nothrow @safe @nogc
if (is(T == int) || is(T == long))
{
    SignedDyadicCoordinate result;

    if (value == 0)
    {
        result.sign = 0;
        return result;
    }

    result.sign =
        value < 0
            ? -1
            : 1;

    const ulong magnitude =
        integralMagnitude(value);

    assert(magnitude != 0);

    setShiftedUnsigned64(
        result.magnitude,
        magnitude,
        dyadicScaleShift
    );

    return result;
}


/**
 * Exact conversion of a supported geo-d scalar coordinate into the
 * common dyadic coordinate representation.
 *
 * All returned values use the same scale:
 *
 *     integer * 2^-1074
 *
 * Supported source types:
 *
 *     int
 *     long
 *     float
 *     double
 *
 * Integral inputs are embedded directly and never pass through a
 * floating-point representation.
 *
 * Every finite binary32 value is exactly representable as binary64, so
 * float conversion may safely reuse the binary64 decoder.
 *
 * Floating inputs must be finite.
 */
SignedDyadicCoordinate decodeDyadicCoordinate(T)(
    T value
)
    pure nothrow @safe @nogc
if (
    is(T == int) ||
    is(T == long) ||
    is(T == float) ||
    is(T == double)
)
{
    static if (
        is(T == int) ||
        is(T == long)
    )
    {
        return decodeIntegralCoordinate(
            value
        );
    }
    else static if (is(T == float))
    {
        assert(isFinite(value));

        /*
         * binary32 -> binary64 is exact.
         */
        return decodeBinary64Coordinate(
            cast(double) value
        );
    }
    else
    {
        assert(isFinite(value));

        return decodeBinary64Coordinate(
            value
        );
    }
}


/**
 * Exact signed difference between two dyadic coordinates.
 *
 * The magnitude uses the same fixed width as the source coordinates.
 * This is sufficient for the complete supported coordinate domain:
 * subtraction of two signed binary64 coordinates may require one bit
 * more than a single coordinate, but still fits in 66 32-bit limbs.
 */
struct SignedDyadicDifference
{
    int sign;
    DyadicCoordinateMagnitude magnitude;
}


/**
 * Computes the exact mathematical difference:
 *
 *     lhs - rhs
 *
 * without signed overflow or floating-point arithmetic.
 *
 * The result remains expressed in units of 2^-1074.
 */
SignedDyadicDifference subtractDyadicCoordinates(
    ref const SignedDyadicCoordinate lhs,
    ref const SignedDyadicCoordinate rhs
)
    pure nothrow @safe @nogc
{
    SignedDyadicDifference result;

    if (lhs.sign == 0)
    {
        result.sign =
            -rhs.sign;

        result.magnitude =
            rhs.magnitude;

        return result;
    }

    if (rhs.sign == 0)
    {
        result.sign =
            lhs.sign;

        result.magnitude =
            lhs.magnitude;

        return result;
    }

    if (lhs.sign != rhs.sign)
    {
        /*
         * Opposite signs turn subtraction into magnitude addition:
         *
         *     (+A) - (-B) = +(A + B)
         *     (-A) - (+B) = -(A + B)
         */
        result.sign =
            lhs.sign;

        result.magnitude =
            addUnsigned(
                lhs.magnitude,
                rhs.magnitude
            );

        return result;
    }

    const int comparison =
        compareUnsigned(
            lhs.magnitude,
            rhs.magnitude
        );

    if (comparison == 0)
    {
        result.sign = 0;
        return result;
    }

    if (comparison > 0)
    {
        result.sign =
            lhs.sign;

        result.magnitude =
            subtractUnsigned(
                lhs.magnitude,
                rhs.magnitude
            );
    }
    else
    {
        result.sign =
            -lhs.sign;

        result.magnitude =
            subtractUnsigned(
                rhs.magnitude,
                lhs.magnitude
            );
    }

    return result;
}


/*
 * A product of two coordinate differences needs at most twice the
 * coordinate-difference width.
 *
 * The numerical scale is:
 *
 *     2^-1074 * 2^-1074
 *   = 2^-2148
 */
enum size_t dyadicProductLimbs =
    2 * dyadicCoordinateLimbs;


alias DyadicProductMagnitude =
    UIntFixed!dyadicProductLimbs;


/**
 * Exact signed product of two dyadic coordinate differences.
 *
 * The magnitude is represented in units of 2^-2148.
 */
struct SignedDyadicProduct
{
    int sign;
    DyadicProductMagnitude magnitude;
}


/**
 * Computes the exact mathematical product:
 *
 *     lhs * rhs
 *
 * without floating-point arithmetic.
 */
SignedDyadicProduct multiplyDyadicDifferences(
    ref const SignedDyadicDifference lhs,
    ref const SignedDyadicDifference rhs
)
    pure nothrow @safe @nogc
{
    SignedDyadicProduct result;

    if (
        lhs.sign == 0 ||
        rhs.sign == 0
    )
    {
        result.sign = 0;
        return result;
    }

    result.sign =
        lhs.sign == rhs.sign
            ? 1
            : -1;

    result.magnitude =
        multiplyUnsigned(
            lhs.magnitude,
            rhs.magnitude
        );

    return result;
}


/**
 * Exact signed difference of two dyadic products.
 *
 * The mathematical result is:
 *
 *     lhs - rhs
 *
 * Products generated by multiplyDyadicDifferences() occupy at most
 * the exact binary64-coordinate product range. Their difference still
 * fits in DyadicProductMagnitude.
 *
 * The result remains expressed in units of 2^-2148.
 */
SignedDyadicProduct subtractDyadicProducts(
    ref const SignedDyadicProduct lhs,
    ref const SignedDyadicProduct rhs
)
    pure nothrow @safe @nogc
{
    SignedDyadicProduct result;

    if (lhs.sign == 0)
    {
        result.sign =
            -rhs.sign;

        result.magnitude =
            rhs.magnitude;

        return result;
    }

    if (rhs.sign == 0)
    {
        result.sign =
            lhs.sign;

        result.magnitude =
            lhs.magnitude;

        return result;
    }

    if (lhs.sign != rhs.sign)
    {
        /*
         * Opposite signs turn subtraction into magnitude addition:
         *
         *     (+A) - (-B) = +(A + B)
         *     (-A) - (+B) = -(A + B)
         */
        result.sign =
            lhs.sign;

        result.magnitude =
            addUnsigned(
                lhs.magnitude,
                rhs.magnitude
            );

        return result;
    }

    const int comparison =
        compareUnsigned(
            lhs.magnitude,
            rhs.magnitude
        );

    if (comparison == 0)
    {
        result.sign = 0;
        return result;
    }

    if (comparison > 0)
    {
        result.sign =
            lhs.sign;

        result.magnitude =
            subtractUnsigned(
                lhs.magnitude,
                rhs.magnitude
            );
    }
    else
    {
        result.sign =
            -lhs.sign;

        result.magnitude =
            subtractUnsigned(
                rhs.magnitude,
                lhs.magnitude
            );
    }

    return result;
}


@safe unittest
{
    const auto zeroCoordinate =
        decodeDyadicCoordinate(0);

    const auto oneCoordinate =
        decodeDyadicCoordinate(1);

    const auto twoCoordinate =
        decodeDyadicCoordinate(2);

    const auto threeCoordinate =
        decodeDyadicCoordinate(3);

    const auto fiveCoordinate =
        decodeDyadicCoordinate(5);

    const auto minusTwoCoordinate =
        decodeDyadicCoordinate(-2);

    const auto oneDifference =
        subtractDyadicCoordinates(
            oneCoordinate,
            zeroCoordinate
        );

    const auto twoDifference =
        subtractDyadicCoordinates(
            twoCoordinate,
            zeroCoordinate
        );

    const auto threeDifference =
        subtractDyadicCoordinates(
            threeCoordinate,
            zeroCoordinate
        );

    const auto fiveDifference =
        subtractDyadicCoordinates(
            fiveCoordinate,
            zeroCoordinate
        );

    const auto minusTwoDifference =
        subtractDyadicCoordinates(
            minusTwoCoordinate,
            zeroCoordinate
        );

    /*
     * 3 - 2 = 1
     */
    {
        const auto three =
            multiplyDyadicDifferences(
                threeDifference,
                oneDifference
            );

        const auto two =
            multiplyDyadicDifferences(
                twoDifference,
                oneDifference
            );

        const auto expected =
            multiplyDyadicDifferences(
                oneDifference,
                oneDifference
            );

        const auto result =
            subtractDyadicProducts(
                three,
                two
            );

        assert(result.sign == 1);

        assert(
            result.magnitude.limb ==
            expected.magnitude.limb
        );
    }


    /*
     * 2 - 3 = -1
     */
    {
        const auto two =
            multiplyDyadicDifferences(
                twoDifference,
                oneDifference
            );

        const auto three =
            multiplyDyadicDifferences(
                threeDifference,
                oneDifference
            );

        const auto expected =
            multiplyDyadicDifferences(
                oneDifference,
                oneDifference
            );

        const auto result =
            subtractDyadicProducts(
                two,
                three
            );

        assert(result.sign == -1);

        assert(
            result.magnitude.limb ==
            expected.magnitude.limb
        );
    }


    /*
     * Equal products produce canonical zero.
     */
    {
        const auto value =
            multiplyDyadicDifferences(
                threeDifference,
                twoDifference
            );

        const auto result =
            subtractDyadicProducts(
                value,
                value
            );

        assert(result.sign == 0);
        assert(result.magnitude.isZero);
    }


    /*
     * 3 - (-2) = 5
     */
    {
        const auto positive =
            multiplyDyadicDifferences(
                threeDifference,
                oneDifference
            );

        const auto negative =
            multiplyDyadicDifferences(
                minusTwoDifference,
                oneDifference
            );

        const auto expected =
            multiplyDyadicDifferences(
                fiveDifference,
                oneDifference
            );

        const auto result =
            subtractDyadicProducts(
                positive,
                negative
            );

        assert(result.sign == 1);

        assert(
            result.magnitude.limb ==
            expected.magnitude.limb
        );
    }


    /*
     * Zero behaves as the additive identity.
     */
    {
        const auto value =
            multiplyDyadicDifferences(
                threeDifference,
                twoDifference
            );

        SignedDyadicProduct zero;

        const auto left =
            subtractDyadicProducts(
                zero,
                value
            );

        assert(left.sign == -1);

        assert(
            left.magnitude.limb ==
            value.magnitude.limb
        );

        const auto right =
            subtractDyadicProducts(
                value,
                zero
            );

        assert(right.sign == 1);

        assert(
            right.magnitude.limb ==
            value.magnitude.limb
        );
    }


    /*
     * The largest binary64 coordinate difference can still participate
     * in an exact product difference without floating-point overflow.
     */
    {
        const auto maximum =
            decodeDyadicCoordinate(
                double.max
            );

        const auto minimum =
            decodeDyadicCoordinate(
                -double.max
            );

        const auto extremeDifference =
            subtractDyadicCoordinates(
                maximum,
                minimum
            );

        const auto extremeProduct =
            multiplyDyadicDifferences(
                extremeDifference,
                extremeDifference
            );

        SignedDyadicProduct zero;

        const auto result =
            subtractDyadicProducts(
                extremeProduct,
                zero
            );

        assert(result.sign == 1);

        assert(
            result.magnitude.limb ==
            extremeProduct.magnitude.limb
        );
    }
}


@safe unittest
{
    /*
     * Unit differences multiply exactly.
     *
     * Each integer coordinate is embedded at bit 1074, therefore
     * 1 * 1 occupies bit 2148 in the product scale.
     */
    {
        const auto zero =
            decodeDyadicCoordinate(0);

        const auto one =
            decodeDyadicCoordinate(1);

        const auto minusOne =
            decodeDyadicCoordinate(-1);

        const auto positiveDifference =
            subtractDyadicCoordinates(
                one,
                zero
            );

        const auto negativeDifference =
            subtractDyadicCoordinates(
                minusOne,
                zero
            );

        const auto positive =
            multiplyDyadicDifferences(
                positiveDifference,
                positiveDifference
            );

        assert(positive.sign == 1);

        enum uint bit =
            2 * dyadicScaleShift;

        assert(
            (
                positive.magnitude.limb[
                    bit / 32
                ] &
                (1U << (bit % 32))
            ) != 0
        );

        const auto negative =
            multiplyDyadicDifferences(
                positiveDifference,
                negativeDifference
            );

        assert(negative.sign == -1);

        assert(
            negative.magnitude.limb ==
            positive.magnitude.limb
        );
    }


    /*
     * Zero difference produces canonical zero product.
     */
    {
        const auto value =
            decodeDyadicCoordinate(
                double.max
            );

        const auto zeroDifference =
            subtractDyadicCoordinates(
                value,
                value
            );

        const auto zero =
            decodeDyadicCoordinate(0.0);

        const auto nonZeroDifference =
            subtractDyadicCoordinates(
                value,
                zero
            );

        const auto product =
            multiplyDyadicDifferences(
                zeroDifference,
                nonZeroDifference
            );

        assert(product.sign == 0);
        assert(product.magnitude.isZero);
    }


    /*
     * Sign multiplication follows ordinary arithmetic.
     */
    {
        const auto zero =
            decodeDyadicCoordinate(0);

        const auto positiveCoordinate =
            decodeDyadicCoordinate(7);

        const auto negativeCoordinate =
            decodeDyadicCoordinate(-11);

        const auto positive =
            subtractDyadicCoordinates(
                positiveCoordinate,
                zero
            );

        const auto negative =
            subtractDyadicCoordinates(
                negativeCoordinate,
                zero
            );

        const auto pp =
            multiplyDyadicDifferences(
                positive,
                positive
            );

        const auto pn =
            multiplyDyadicDifferences(
                positive,
                negative
            );

        const auto nn =
            multiplyDyadicDifferences(
                negative,
                negative
            );

        assert(pp.sign == 1);
        assert(pn.sign == -1);
        assert(nn.sign == 1);
    }


    /*
     * Extreme binary64 differences remain inside the fixed product
     * width.
     */
    {
        const auto positive =
            decodeDyadicCoordinate(
                double.max
            );

        const auto negative =
            decodeDyadicCoordinate(
                -double.max
            );

        const auto difference =
            subtractDyadicCoordinates(
                positive,
                negative
            );

        const auto product =
            multiplyDyadicDifferences(
                difference,
                difference
            );

        assert(product.sign == 1);
        assert(!product.magnitude.isZero);

        static assert(
            dyadicProductLimbs ==
            132
        );
    }


    /*
     * Ordinary signed subtraction.
     */
    {
        const auto one =
            decodeDyadicCoordinate(1);

        const auto negativeOne =
            decodeDyadicCoordinate(-1);

        const auto expectedTwo =
            decodeDyadicCoordinate(2);

        const auto positive =
            subtractDyadicCoordinates(
                one,
                negativeOne
            );

        assert(positive.sign == 1);

        assert(
            positive.magnitude.limb ==
            expectedTwo.magnitude.limb
        );

        const auto negative =
            subtractDyadicCoordinates(
                negativeOne,
                one
            );

        assert(negative.sign == -1);

        assert(
            negative.magnitude.limb ==
            expectedTwo.magnitude.limb
        );
    }


    /*
     * Equal values produce canonical exact zero.
     */
    {
        const auto value =
            decodeDyadicCoordinate(
                long.max
            );

        const auto difference =
            subtractDyadicCoordinates(
                value,
                value
            );

        assert(difference.sign == 0);
        assert(difference.magnitude.isZero);
    }


    /*
     * Signed zero participates normally.
     */
    {
        const auto zero =
            decodeDyadicCoordinate(0);

        const auto value =
            decodeDyadicCoordinate(17);

        const auto left =
            subtractDyadicCoordinates(
                zero,
                value
            );

        assert(left.sign == -1);

        assert(
            left.magnitude.limb ==
            value.magnitude.limb
        );

        const auto right =
            subtractDyadicCoordinates(
                value,
                zero
            );

        assert(right.sign == 1);

        assert(
            right.magnitude.limb ==
            value.magnitude.limb
        );
    }


    /*
     * Difference may exceed the range of the original integral scalar.
     *
     * This exercises long.min / long.max without signed overflow.
     */
    {
        const auto minimum =
            decodeDyadicCoordinate(
                long.min
            );

        const auto maximum =
            decodeDyadicCoordinate(
                long.max
            );

        const auto difference =
            subtractDyadicCoordinates(
                maximum,
                minimum
            );

        assert(difference.sign == 1);

        /*
         * long.max - long.min = 2^64 - 1.
         *
         * In the common dyadic scale, bits 1074 through 1137 are all
         * set.
         */
        foreach (
            bit;
            dyadicScaleShift ..
            dyadicScaleShift + 64
        )
        {
            assert(
                (
                    difference.magnitude.limb[
                        bit / 32
                    ] &
                    (1U << (bit % 32))
                ) != 0
            );
        }
    }


    /*
     * Difference of opposite maximum binary64 coordinates needs one
     * more significant bit than a single coordinate and still fits the
     * shared fixed-width representation.
     */
    {
        const auto positive =
            decodeDyadicCoordinate(
                double.max
            );

        const auto negative =
            decodeDyadicCoordinate(
                -double.max
            );

        const auto difference =
            subtractDyadicCoordinates(
                positive,
                negative
            );

        assert(difference.sign == 1);
        assert(!difference.magnitude.isZero);

        /*
         * double.max has highest scaled-integer bit 2097.
         * Doubling it produces a result whose highest set bit is 2098.
         */
        enum uint highestBit = 2098;

        assert(
            (
                difference.magnitude.limb[
                    highestBit / 32
                ] &
                (1U << (highestBit % 32))
            ) != 0
        );
    }


    /*
     * Subnormal differences remain exact.
     */
    {
        enum double tiny =
            0x0.0000000000001p-1022;

        enum double twiceTiny =
            0x0.0000000000002p-1022;

        const auto one =
            decodeDyadicCoordinate(
                tiny
            );

        const auto two =
            decodeDyadicCoordinate(
                twiceTiny
            );

        const auto difference =
            subtractDyadicCoordinates(
                two,
                one
            );

        assert(difference.sign == 1);

        assert(
            difference.magnitude.limb ==
            one.magnitude.limb
        );
    }


    /*
     * Coordinates originating from different supported scalar types
     * share exactly the same dyadic number line.
     */
    {
        const auto integer =
            decodeDyadicCoordinate(1L);

        const auto floating =
            decodeDyadicCoordinate(1.0);

        const auto difference =
            subtractDyadicCoordinates(
                integer,
                floating
            );

        assert(difference.sign == 0);
        assert(difference.magnitude.isZero);
    }


    /*
     * Integral unit coordinates occupy bit 1074 in the common scale.
     */
    {
        const auto positive =
            decodeDyadicCoordinate(1);

        const auto negative =
            decodeDyadicCoordinate(-1);

        assert(positive.sign == 1);
        assert(negative.sign == -1);

        enum size_t limbIndex =
            dyadicScaleShift / 32;

        enum uint bitIndex =
            dyadicScaleShift % 32;

        assert(
            positive.magnitude
                .limb[limbIndex] ==
            (1U << bitIndex)
        );

        assert(
            positive.magnitude.limb ==
            negative.magnitude.limb
        );
    }


    /*
     * int.min is converted without signed overflow.
     *
     *     |int.min| = 2^31
     */
    {
        const auto value =
            decodeDyadicCoordinate(
                int.min
            );

        assert(value.sign == -1);

        enum uint bit =
            dyadicScaleShift + 31;

        enum size_t limbIndex =
            bit / 32;

        enum uint bitIndex =
            bit % 32;

        assert(
            value.magnitude
                .limb[limbIndex] ==
            (1U << bitIndex)
        );
    }


    /*
     * long.min is converted without ever evaluating -long.min.
     *
     *     |long.min| = 2^63
     */
    {
        const auto value =
            decodeDyadicCoordinate(
                long.min
            );

        assert(value.sign == -1);

        enum uint bit =
            dyadicScaleShift + 63;

        enum size_t limbIndex =
            bit / 32;

        enum uint bitIndex =
            bit % 32;

        assert(
            value.magnitude
                .limb[limbIndex] ==
            (1U << bitIndex)
        );
    }


    /*
     * long.max remains exact as well.
     */
    {
        const auto value =
            decodeDyadicCoordinate(
                long.max
            );

        assert(value.sign == 1);

        /*
         * long.max contains bits 0 .. 62. After embedding, bit 1074 is
         * the least significant set bit and bit 1136 the greatest.
         */
        enum uint lowest =
            dyadicScaleShift;

        enum uint highest =
            dyadicScaleShift + 62;

        assert(
            (
                value.magnitude.limb[
                    lowest / 32
                ] &
                (1U << (lowest % 32))
            ) != 0
        );

        assert(
            (
                value.magnitude.limb[
                    highest / 32
                ] &
                (1U << (highest % 32))
            ) != 0
        );
    }


    /*
     * Integral zero normalizes to the same signed-zero representation
     * as floating zero.
     */
    {
        const auto integerZero =
            decodeDyadicCoordinate(0);

        const auto longZero =
            decodeDyadicCoordinate(0L);

        const auto doubleZero =
            decodeDyadicCoordinate(0.0);

        assert(integerZero.sign == 0);
        assert(longZero.sign == 0);
        assert(doubleZero.sign == 0);

        assert(
            integerZero.magnitude.isZero
        );

        assert(
            longZero.magnitude.isZero
        );

        assert(
            doubleZero.magnitude.isZero
        );
    }


    /*
     * binary32 conversion is exactly equivalent to exact promotion to
     * binary64 followed by binary64 decoding.
     */
    {
        enum float[] values = [
            1.0f,
            -1.0f,
            0x1.000002p+0f,
            float.max,
            -float.max,
            0x1p-149f
        ];

        static foreach (value; values)
        {{
            const auto fromFloat =
                decodeDyadicCoordinate(
                    value
                );

            const auto fromDouble =
                decodeBinary64Coordinate(
                    cast(double) value
                );

            assert(
                fromFloat.sign ==
                fromDouble.sign
            );

            assert(
                fromFloat.magnitude.limb ==
                fromDouble.magnitude.limb
            );
        }}
    }


    /*
     * The smallest positive binary32 subnormal corresponds exactly to:
     *
     *     2^-149
     *
     * In the 2^-1074 common scale this is bit:
     *
     *     1074 - 149 = 925
     */
    {
        enum float smallest =
            0x1p-149f;

        const auto value =
            decodeDyadicCoordinate(
                smallest
            );

        assert(value.sign == 1);

        enum uint bit =
            dyadicScaleShift - 149;

        assert(
            value.magnitude.limb[
                bit / 32
            ] ==
            (1U << (bit % 32))
        );
    }


    /*
     * Generic binary64 conversion is exactly the existing decoder.
     */
    {
        enum double[] values = [
            1.0,
            -1.0,
            0x1.0000000000001p+0,
            double.max,
            -double.max,
            0x0.0000000000001p-1022
        ];

        static foreach (value; values)
        {{
            const auto generic =
                decodeDyadicCoordinate(
                    value
                );

            const auto direct =
                decodeBinary64Coordinate(
                    value
                );

            assert(
                generic.sign ==
                direct.sign
            );

            assert(
                generic.magnitude.limb ==
                direct.magnitude.limb
            );
        }}
    }


    /*
     * The common converter intentionally excludes real until a
     * platform-aware backend is designed.
     */
    static assert(
        !__traits(
            compiles,
            decodeDyadicCoordinate(
                cast(real) 1
            )
        )
    );


    /*
     * Both signed zeros normalize to exact dyadic zero.
     */
    {
        const auto positive =
            decodeBinary64Coordinate(
                0.0
            );

        const auto negative =
            decodeBinary64Coordinate(
                -0.0
            );

        assert(positive.sign == 0);
        assert(negative.sign == 0);

        assert(
            positive.magnitude.isZero
        );

        assert(
            negative.magnitude.isZero
        );
    }


    /*
     * The smallest positive binary64 subnormal is exactly one unit of
     * the common 2^-1074 coordinate scale.
     */
    {
        enum double smallest =
            0x0.0000000000001p-1022;

        const auto value =
            decodeBinary64Coordinate(
                smallest
            );

        assert(value.sign == 1);
        assert(
            value.magnitude.limb[0] ==
            1
        );

        foreach (
            index;
            1 .. dyadicCoordinateLimbs
        )
        {
            assert(
                value.magnitude
                    .limb[index] == 0
            );
        }
    }


    /*
     * Sign is separated from magnitude.
     */
    {
        const auto positive =
            decodeBinary64Coordinate(
                1.0
            );

        const auto negative =
            decodeBinary64Coordinate(
                -1.0
            );

        assert(positive.sign == 1);
        assert(negative.sign == -1);

        assert(
            positive.magnitude.limb ==
            negative.magnitude.limb
        );
    }


    /*
     * The complete finite binary64 range fits in the fixed coordinate
     * width.
     */
    {
        const auto maximum =
            decodeBinary64Coordinate(
                double.max
            );

        const auto minimum =
            decodeBinary64Coordinate(
                -double.max
            );

        assert(maximum.sign == 1);
        assert(minimum.sign == -1);

        assert(
            !maximum.magnitude.isZero
        );

        assert(
            maximum.magnitude.limb ==
            minimum.magnitude.limb
        );
    }
}
