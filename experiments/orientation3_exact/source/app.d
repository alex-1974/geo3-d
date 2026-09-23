module app;

import geo3 : Point3;

import std.bigint : BigInt;
import std.stdio : writeln;


/*
 * Research-only fixed-width unsigned integer.
 *
 * Little-endian base-2^32 limbs.
 */
struct UIntFixed(size_t Limbs)
if (Limbs > 0)
{
    uint[Limbs] limb;
}


int compareUnsigned(size_t Limbs)(
    ref const UIntFixed!Limbs lhs,
    ref const UIntFixed!Limbs rhs
)
    pure nothrow @safe @nogc
{
    for (size_t i = Limbs; i != 0; --i)
    {
        const size_t index = i - 1;

        if (lhs.limb[index] < rhs.limb[index])
            return -1;

        if (lhs.limb[index] > rhs.limb[index])
            return 1;
    }

    return 0;
}


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

    foreach (i; 0 .. ValueLimbs)
    {
        const ulong sum =
            cast(ulong) target.limb[i] +
            cast(ulong) value.limb[i] +
            carry;

        target.limb[i] = cast(uint) sum;
        carry = sum >> 32;
    }

    size_t index = ValueLimbs;

    while (carry != 0 && index < TargetLimbs)
    {
        const ulong sum =
            cast(ulong) target.limb[index] +
            carry;

        target.limb[index] = cast(uint) sum;
        carry = sum >> 32;

        ++index;
    }

    return carry == 0;
}


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
            const size_t index = i + j;

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
            assert(index < LhsLimbs + RhsLimbs);

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


struct SignedDiff(size_t Limbs)
{
    int sign;
    UIntFixed!Limbs magnitude;
}


uint unsignedMagnitude(int value)
    pure nothrow @safe @nogc
{
    if (value >= 0)
        return cast(uint) value;

    return cast(uint)(-(value + 1)) + 1U;
}


ulong unsignedMagnitude(long value)
    pure nothrow @safe @nogc
{
    if (value >= 0)
        return cast(ulong) value;

    return cast(ulong)(-(value + 1)) + 1UL;
}


SignedDiff!1 signedDifference(int lhs, int rhs)
    pure nothrow @safe @nogc
{
    SignedDiff!1 result;

    if (lhs == rhs)
        return result;

    result.sign =
        lhs > rhs ? 1 : -1;

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
            cast(uint)(lhs - rhs);
    }
    else
    {
        magnitude =
            cast(uint)(rhs - lhs);
    }

    result.magnitude.limb[0] =
        magnitude;

    return result;
}


SignedDiff!2 signedDifference(long lhs, long rhs)
    pure nothrow @safe @nogc
{
    SignedDiff!2 result;

    if (lhs == rhs)
        return result;

    result.sign =
        lhs > rhs ? 1 : -1;

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
            cast(ulong)(lhs - rhs);
    }
    else
    {
        magnitude =
            cast(ulong)(rhs - lhs);
    }

    result.magnitude.limb[0] =
        cast(uint) magnitude;

    result.magnitude.limb[1] =
        cast(uint)(magnitude >> 32);

    return result;
}


void accumulateTerm(
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
if (BucketLimbs >= DiffLimbs * 3)
{
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

    UIntFixed!BucketLimbs widened;

    foreach (i; 0 .. xyz.limb.length)
        widened.limb[i] = xyz.limb[i];

    bool success;

    if (sign > 0)
        success = addAssign(positive, widened);
    else
        success = addAssign(negative, widened);

    assert(success);
}


int orientationExactImpl(
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
       is(T == int)
    || is(T == long)
)
{
    const auto ux = signedDifference(b.x, a.x);
    const auto uy = signedDifference(b.y, a.y);
    const auto uz = signedDifference(b.z, a.z);

    const auto vx = signedDifference(c.x, a.x);
    const auto vy = signedDifference(c.y, a.y);
    const auto vz = signedDifference(c.z, a.z);

    const auto wx = signedDifference(d.x, a.x);
    const auto wy = signedDifference(d.y, a.y);
    const auto wz = signedDifference(d.z, a.z);

    UIntFixed!BucketLimbs positive;
    UIntFixed!BucketLimbs negative;

    accumulateTerm(
        positive, negative,
        +1,
        ux, vy, wz
    );

    accumulateTerm(
        positive, negative,
        +1,
        uy, vz, wx
    );

    accumulateTerm(
        positive, negative,
        +1,
        uz, vx, wy
    );

    accumulateTerm(
        positive, negative,
        -1,
        uz, vy, wx
    );

    accumulateTerm(
        positive, negative,
        -1,
        uy, vx, wz
    );

    accumulateTerm(
        positive, negative,
        -1,
        ux, vz, wy
    );

    return compareUnsigned(
        positive,
        negative
    );
}


int orientationExact(
    Point3!int a,
    Point3!int b,
    Point3!int c,
    Point3!int d
)
    pure nothrow @safe @nogc
{
    return orientationExactImpl!(
        int,
        1,
        4
    )(a, b, c, d);
}


int orientationExact(
    Point3!long a,
    Point3!long b,
    Point3!long c,
    Point3!long d
)
    pure nothrow @safe @nogc
{
    return orientationExactImpl!(
        long,
        2,
        7
    )(a, b, c, d);
}


/*
 * Independent arbitrary-precision oracle.
 */
int oracleOrientation(T)(
    Point3!T a,
    Point3!T b,
    Point3!T c,
    Point3!T d
)
    @safe
if (
       is(T == int)
    || is(T == long)
)
{
    const BigInt ux =
        BigInt(b.x) - BigInt(a.x);

    const BigInt uy =
        BigInt(b.y) - BigInt(a.y);

    const BigInt uz =
        BigInt(b.z) - BigInt(a.z);

    const BigInt vx =
        BigInt(c.x) - BigInt(a.x);

    const BigInt vy =
        BigInt(c.y) - BigInt(a.y);

    const BigInt vz =
        BigInt(c.z) - BigInt(a.z);

    const BigInt wx =
        BigInt(d.x) - BigInt(a.x);

    const BigInt wy =
        BigInt(d.y) - BigInt(a.y);

    const BigInt wz =
        BigInt(d.z) - BigInt(a.z);

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


void checkExpected(T)(
    int expected,
    Point3!T a,
    Point3!T b,
    Point3!T c,
    Point3!T d
)
    @safe
if (
       is(T == int)
    || is(T == long)
)
{
    assert(expected >= -1);
    assert(expected <= 1);

    assert(
        oracleOrientation(a, b, c, d) ==
        expected
    );

    assert(
        orientationExact(a, b, c, d) ==
        expected
    );

    checkCase(a, b, c, d);
}


void checkCase(T)(
    Point3!T a,
    Point3!T b,
    Point3!T c,
    Point3!T d
)
    @safe
if (
       is(T == int)
    || is(T == long)
)
{
    const int expected =
        oracleOrientation(a, b, c, d);

    const int actual =
        orientationExact(a, b, c, d);

    assert(actual == expected);

    /*
     * One transposition reverses sign.
     */
    assert(
        orientationExact(a, c, b, d) ==
        -actual
    );

    /*
     * A four-cycle is odd and reverses sign.
     */
    assert(
        orientationExact(b, c, d, a) ==
        -actual
    );

    /*
     * Two transpositions preserve sign.
     */
    assert(
        orientationExact(b, a, d, c) ==
        actual
    );
}


ulong nextRandom(ref ulong state)
    pure nothrow @safe @nogc
{
    assert(state != 0);

    state ^= state << 13;
    state ^= state >> 7;
    state ^= state << 17;

    return state;
}


Point3!int randomIntPoint(ref ulong state)
    pure nothrow @safe @nogc
{
    return Point3!int(
        cast(int) nextRandom(state),
        cast(int) nextRandom(state),
        cast(int) nextRandom(state)
    );
}


Point3!long randomLongPoint(ref ulong state)
    pure nothrow @safe @nogc
{
    return Point3!long(
        cast(long) nextRandom(state),
        cast(long) nextRandom(state),
        cast(long) nextRandom(state)
    );
}


int signFromMask(
    size_t mask,
    size_t bit
)
    pure nothrow @safe @nogc
{
    return (
        mask &
        (cast(size_t) 1 << bit)
    ) != 0
        ? 1
        : -1;
}


void verifyLeibnizSignBuckets()
    pure nothrow @safe @nogc
{
    size_t splitFiveOne = 0;
    size_t splitThreeThree = 0;
    size_t splitOneFive = 0;

    foreach (mask; 0 .. 512)
    {
        int[9] sign;

        foreach (i; 0 .. sign.length)
        {
            sign[i] =
                signFromMask(mask, i);
        }

        const int ux = sign[0];
        const int uy = sign[1];
        const int uz = sign[2];

        const int vx = sign[3];
        const int vy = sign[4];
        const int vz = sign[5];

        const int wx = sign[6];
        const int wy = sign[7];
        const int wz = sign[8];

        int[6] termSign;

        termSign[0] =
            ux * vy * wz;

        termSign[1] =
            uy * vz * wx;

        termSign[2] =
            uz * vx * wy;

        termSign[3] =
            -uz * vy * wx;

        termSign[4] =
            -uy * vx * wz;

        termSign[5] =
            -ux * vz * wy;

        size_t positive = 0;
        size_t negative = 0;
        int signProduct = 1;

        foreach (value; termSign)
        {
            assert(
                value == -1 ||
                value == 1
            );

            signProduct *= value;

            if (value > 0)
                ++positive;
            else
                ++negative;
        }

        assert(
            positive + negative == 6
        );

        /*
         * Product of all six Leibniz term signs is always negative.
         * Therefore the number of negative terms is always odd.
         */
        assert(signProduct == -1);

        if (
            positive == 5 &&
            negative == 1
        )
        {
            ++splitFiveOne;
        }
        else if (
            positive == 3 &&
            negative == 3
        )
        {
            ++splitThreeThree;
        }
        else if (
            positive == 1 &&
            negative == 5
        )
        {
            ++splitOneFive;
        }
        else
        {
            assert(0);
        }
    }

    /*
     * Exhaustive distribution over all 2^9 non-zero component-sign
     * assignments.
     */
    assert(splitFiveOne == 96);
    assert(splitThreeThree == 320);
    assert(splitOneFive == 96);
}


void main()
    @safe
{
    alias PI = Point3!int;
    alias PL = Point3!long;

    verifyLeibnizSignBuckets();

    /*
     * Canonical public sign convention:
     *
     * det(b-a, c-a, d-a) > 0.
     */
    checkExpected(
        +1,
        PI(0, 0, 0),
        PI(1, 0, 0),
        PI(0, 1, 0),
        PI(0, 0, 1)
    );

    /*
     * Exactly coplanar.
     */
    checkExpected(
        0,
        PI(0, 0, 0),
        PI(1, 0, 0),
        PI(0, 1, 0),
        PI(1, 1, 0)
    );

    /*
     * Degenerate / repeated points.
     */
    checkExpected(
        0,
        PI(7, -2, 9),
        PI(7, -2, 9),
        PI(3, 5, 11),
        PI(-8, 4, 6)
    );

    /*
     * Complete int span.
     */
    checkExpected(
        +1,
        PI(int.min, int.min, int.min),
        PI(int.max, int.min, int.min),
        PI(int.min, int.max, int.min),
        PI(int.min, int.min, int.max)
    );

    /*
     * Complete long span.
     *
     * The determinant magnitude is approximately (2^64)^3 and cannot
     * be represented by native 64- or 128-bit signed arithmetic.
     */
    checkExpected(
        +1,
        PL(long.min, long.min, long.min),
        PL(long.max, long.min, long.min),
        PL(long.min, long.max, long.min),
        PL(long.min, long.min, long.max)
    );

    /*
     * Exact coplanarity across complete X/Y coordinate spans.
     *
     * All four points lie in z = long.min.
     */
    checkExpected(
        0,
        PL(long.min, long.min, long.min),
        PL(long.max, long.min, long.min),
        PL(long.min, long.max, long.min),
        PL(long.max, long.max, long.min)
    );

    /*
     * Large cancellation-heavy coplanarity with no zero direction
     * components.
     *
     * d-a == (b-a) + (c-a), so the three difference vectors are
     * linearly dependent. Individual triple products are far wider
     * than native 128-bit arithmetic.
     */
    enum long m =
        (1L << 61) - 1;

    checkExpected(
        0,
        PL(0, 0, 0),
        PL(m, 2 * m, 3 * m),
        PL(2 * m, -m, m),
        PL(3 * m, m, 4 * m)
    );

    ulong state =
        0x6a09_e667_f3bc_c909UL;

    enum randomCases = 50_000;

    foreach (_; 0 .. randomCases)
    {
        checkCase(
            randomIntPoint(state),
            randomIntPoint(state),
            randomIntPoint(state),
            randomIntPoint(state)
        );
    }

    foreach (_; 0 .. randomCases)
    {
        checkCase(
            randomLongPoint(state),
            randomLongPoint(state),
            randomLongPoint(state),
            randomLongPoint(state)
        );
    }

    writeln(
        "orientation3 exact prototype: PASS (",
        randomCases,
        " int + ",
        randomCases,
        " long random cases)"
    );
}
