/**
 * Exact fixed-width dyadic representation of finite binary64 coordinates.
 *
 * INTERNAL IMPLEMENTATION MODULE.
 *
 * Every finite IEEE binary64 value is represented exactly as:
 *
 *     integer * 2^-1074
 *
 * 66 32-bit limbs provide 2112 bits. This covers every finite binary64
 * coordinate and every exact difference between two such coordinates.
 *
 * This module contains only the functionality required by Orientation3.
 *
 * Provenance:
 *     Derived from geo-d v2.0.0:
 *
 *         repository: https://github.com/alex-1974/geo-d
 *         tag:        v2.0.0
 *         commit:     83c974e0ca018bb378ce65a08ef9d3178decf9ad
 *         source:     source/geo/internal/dyadic.d
 *         blob:       554a657f271da803dd6e524117d1e5bc6986090a
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
module geo3.internal.dyadic;

import geo3.internal.fixed_uint :
    UIntFixed,
    addAssign,
    compareUnsigned,
    multiplyUnsigned,
    subtractUnsigned;

import std.bitmanip :
    DoubleRep;

import std.math.traits :
    isFinite;


/**
 * Number of 32-bit limbs required by every finite binary64 coordinate or
 * exact coordinate difference in the common 2^-1074 scale.
 */
enum size_t dyadicCoordinateLimbs =
    66;


/*
 * Every magnitude is represented as an integer multiple of 2^-1074.
 */
private enum uint dyadicScaleShift =
    1074;


alias DyadicCoordinateMagnitude =
    UIntFixed!dyadicCoordinateLimbs;


/**
 * Exact signed binary64 coordinate.
 */
struct SignedDyadicCoordinate
{
    int sign;
    DyadicCoordinateMagnitude magnitude;
}


/**
 * Exact signed difference between two binary64 coordinates.
 */
struct SignedDyadicDifference
{
    int sign;
    DyadicCoordinateMagnitude magnitude;
}


/**
 * Number of limbs required by a product of two coordinate differences.
 */
enum size_t dyadicProductLimbs =
    2 * dyadicCoordinateLimbs;


alias DyadicProductMagnitude =
    UIntFixed!dyadicProductLimbs;


/**
 * Exact signed product of two coordinate differences.
 */
struct SignedDyadicProduct
{
    int sign;
    DyadicProductMagnitude magnitude;
}


/*
 * Places a non-zero 64-bit mantissa at the requested bit position.
 *
 * The destination must initially be zero.
 */
private void setShiftedUnsigned64(
    ref DyadicCoordinateMagnitude result,
    ulong mantissa,
    uint shift
)
    pure nothrow @safe @nogc
{
    assert(mantissa != 0);

    /*
     * Maximum normal binary64 shift:
     *
     *     rawExponent - 1
     *     = 2046 - 1
     *     = 2045
     */
    assert(shift <= 2045);


    const size_t base =
        shift / 32;

    const uint offset =
        shift % 32;

    const uint lower =
        cast(uint) mantissa;

    const uint upper =
        cast(uint)(
            mantissa >> 32
        );


    const ulong shiftedLower =
        cast(ulong) lower <<
        offset;

    result.limb[base] |=
        cast(uint) shiftedLower;


    if (
        base + 1 <
        dyadicCoordinateLimbs
    )
    {
        result.limb[base + 1] |=
            cast(uint)(
                shiftedLower >> 32
            );
    }


    if (upper != 0)
    {
        const ulong shiftedUpper =
            cast(ulong) upper <<
            offset;

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
                (
                    shiftedUpper >>
                    32
                ) == 0
            );
        }
    }
}


/**
 * Decodes one finite binary64 coordinate exactly.
 *
 * Subnormals:
 *
 *     value = fraction * 2^-1074
 *
 * Normals:
 *
 *     value =
 *         (2^52 + fraction)
 *         * 2^(rawExponent - 1075)
 *
 * Therefore the exact integer in the common 2^-1074 scale is:
 *
 *     mantissa << (rawExponent - 1)
 *
 * Positive and negative zero both normalize to sign == 0.
 */
SignedDyadicCoordinate decodeBinary64Coordinate(
    double value
)
    pure nothrow @safe @nogc
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


    SignedDyadicCoordinate result;


    if (mantissa == 0)
    {
        result.sign =
            0;

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


/**
 * Exact mathematical difference:
 *
 *     lhs - rhs
 *
 * The result remains in the common 2^-1074 scale.
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
         * (+A) - (-B) = +(A+B)
         * (-A) - (+B) = -(A+B)
         */
        result.sign =
            lhs.sign;

        result.magnitude =
            lhs.magnitude;

        const bool success =
            addAssign(
                result.magnitude,
                rhs.magnitude
            );

        /*
         * The exact difference of two finite binary64 coordinates requires
         * at most bit 2098 and therefore fits in 66 limbs.
         */
        assert(success);

        return result;
    }


    const int comparison =
        compareUnsigned(
            lhs.magnitude,
            rhs.magnitude
        );


    if (comparison == 0)
    {
        result.sign =
            0;

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


/**
 * Exact product of two dyadic coordinate differences.
 *
 * The result is represented in units of:
 *
 *     2^-2148
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
        result.sign =
            0;

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


static assert(
    dyadicCoordinateLimbs ==
    66
);

static assert(
    dyadicCoordinateLimbs * 32 ==
    2112
);

static assert(
    dyadicProductLimbs ==
    132
);


@safe unittest
{
    enum double minSubnormal =
        0x0.0000000000001p-1022;


    /*
     * Smallest binary64 subnormal is exactly one unit of the common scale.
     */
    const auto smallest =
        decodeBinary64Coordinate(
            minSubnormal
        );

    assert(
        smallest.sign ==
        1
    );

    assert(
        smallest.magnitude.limb[0] ==
        1
    );

    foreach (
        index;
        1 .. dyadicCoordinateLimbs
    )
    {
        assert(
            smallest.magnitude
                .limb[index] == 0
        );
    }


    /*
     * Signed zeros normalize to the same mathematical zero.
     */
    const auto positiveZero =
        decodeBinary64Coordinate(
            0.0
        );

    const auto negativeZero =
        decodeBinary64Coordinate(
            -0.0
        );

    assert(
        positiveZero.sign ==
        0
    );

    assert(
        negativeZero.sign ==
        0
    );


    /*
     * Opposite maximum finite coordinates require one additional bit in
     * their exact difference and still fit in 66 limbs.
     */
    const auto maximum =
        decodeBinary64Coordinate(
            double.max
        );

    const auto minimum =
        decodeBinary64Coordinate(
            -double.max
        );

    const auto completeSpan =
        subtractDyadicCoordinates(
            maximum,
            minimum
        );

    assert(
        completeSpan.sign ==
        1
    );

    enum uint highestDifferenceBit =
        2098;

    assert(
        (
            completeSpan.magnitude.limb[
                highestDifferenceBit /
                32
            ] &
            (
                1U <<
                (
                    highestDifferenceBit %
                    32
                )
            )
        ) != 0
    );


    /*
     * Subnormal subtraction remains exact.
     */
    enum double twiceMinSubnormal =
        0x0.0000000000002p-1022;

    const auto twiceSmallest =
        decodeBinary64Coordinate(
            twiceMinSubnormal
        );

    const auto subnormalDifference =
        subtractDyadicCoordinates(
            twiceSmallest,
            smallest
        );

    assert(
        subnormalDifference.sign ==
        1
    );

    assert(
        subnormalDifference.magnitude
            .limb[0] == 1
    );


    /*
     * Unit-coordinate differences produce their product at bit 2148.
     */
    const auto zero =
        decodeBinary64Coordinate(
            0.0
        );

    const auto one =
        decodeBinary64Coordinate(
            1.0
        );

    const auto oneDifference =
        subtractDyadicCoordinates(
            one,
            zero
        );

    const auto unitProduct =
        multiplyDyadicDifferences(
            oneDifference,
            oneDifference
        );

    assert(
        unitProduct.sign ==
        1
    );

    enum uint productBit =
        2 * dyadicScaleShift;

    assert(
        (
            unitProduct.magnitude.limb[
                productBit / 32
            ] &
            (
                1U <<
                (
                    productBit %
                    32
                )
            )
        ) != 0
    );
}
