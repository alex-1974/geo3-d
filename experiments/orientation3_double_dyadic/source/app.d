module app;

import geo3 : Point3;

import geo.internal.dyadic :
    SignedDyadicDifference,
    decodeBinary64Coordinate,
    dyadicCoordinateLimbs,
    multiplyDyadicDifferences,
    subtractDyadicCoordinates;

import geo.internal.fixed_uint :
    UIntFixed,
    addUnsigned,
    compareUnsigned,
    multiplyUnsigned;

import std.bigint : BigInt;
import std.bitmanip : DoubleRep;
import std.math.traits : isFinite;
import std.stdio : writeln;


alias P = Point3!double;


enum size_t tripleProductLimbs =
    3 * dyadicCoordinateLimbs;

static assert(dyadicCoordinateLimbs == 66);
static assert(tripleProductLimbs == 198);
static assert(tripleProductLimbs * 32 == 6336);

alias TripleMagnitude =
    UIntFixed!tripleProductLimbs;


/*
 * Accumulate one exact signed triple product:
 *
 *     coefficient * x * y * z
 *
 * into positive or negative unsigned determinant buckets.
 */
void accumulateTriple(
    ref TripleMagnitude positive,
    ref TripleMagnitude negative,
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

    assert(xy.sign != 0);

    const TripleMagnitude xyz =
        multiplyUnsigned(
            xy.magnitude,
            z.magnitude
        );

    const int sign =
        coefficient *
        xy.sign *
        z.sign;

    if (sign > 0)
    {
        positive =
            addUnsigned(
                positive,
                xyz
            );
    }
    else
    {
        negative =
            addUnsigned(
                negative,
                xyz
            );
    }
}


/*
 * Full-range exact orient3d for every finite binary64 input.
 *
 * Computes:
 *
 *     sign(det(b-a, c-a, d-a))
 *
 * No floating-point arithmetic occurs after coordinate decoding.
 */
int orientationDyadicExact3(
    P a,
    P b,
    P c,
    P d
)
    pure nothrow @safe @nogc
{
    assert(a.isFinite);
    assert(b.isFinite);
    assert(c.isFinite);
    assert(d.isFinite);


    const auto ax =
        decodeBinary64Coordinate(a.x);

    const auto ay =
        decodeBinary64Coordinate(a.y);

    const auto az =
        decodeBinary64Coordinate(a.z);

    const auto bx =
        decodeBinary64Coordinate(b.x);

    const auto by =
        decodeBinary64Coordinate(b.y);

    const auto bz =
        decodeBinary64Coordinate(b.z);

    const auto cx =
        decodeBinary64Coordinate(c.x);

    const auto cy =
        decodeBinary64Coordinate(c.y);

    const auto cz =
        decodeBinary64Coordinate(c.z);

    const auto dx =
        decodeBinary64Coordinate(d.x);

    const auto dy =
        decodeBinary64Coordinate(d.y);

    const auto dz =
        decodeBinary64Coordinate(d.z);


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


    TripleMagnitude positive;
    TripleMagnitude negative;


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


/*
 * Independent exact binary64 -> integer * 2^-1074 oracle.
 */
BigInt oracleBinary64Integer(double value)
    @safe
{
    assert(isFinite(value));

    DoubleRep representation;
    representation.value = value;

    const ulong fraction =
        representation.fraction;

    const uint rawExponent =
        representation.exponent;

    ulong mantissa;
    uint shift;

    if (rawExponent == 0)
    {
        mantissa = fraction;
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
 * Independent BigInt orient3d oracle.
 */
int oracleOrientation(
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


    const BigInt ux = bx - ax;
    const BigInt uy = by - ay;
    const BigInt uz = bz - az;

    const BigInt vx = cx - ax;
    const BigInt vy = cy - ay;
    const BigInt vz = cz - az;

    const BigInt wx = dx - ax;
    const BigInt wy = dy - ay;
    const BigInt wz = dz - az;


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


void checkCase(
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
        orientationDyadicExact3(
            a,
            b,
            c,
            d
        );

    assert(actual == expected);


    /*
     * One transposition reverses sign.
     */
    assert(
        orientationDyadicExact3(
            a,
            c,
            b,
            d
        ) == -actual
    );


    /*
     * A four-cycle is odd and reverses sign.
     */
    assert(
        orientationDyadicExact3(
            b,
            c,
            d,
            a
        ) == -actual
    );
}


void checkExpected(
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
        ) == expected
    );

    assert(
        orientationDyadicExact3(
            a,
            b,
            c,
            d
        ) == expected
    );

    checkCase(
        a,
        b,
        c,
        d
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


double randomFiniteDouble(ref ulong state)
    pure nothrow @safe @nogc
{
    const ulong bits =
        nextRandom(state);

    DoubleRep representation;
    representation.value = 0.0;

    representation.fraction =
        bits &
        ((1UL << 52) - 1);

    ushort exponent =
        cast(ushort)(
            (bits >> 52) &
            0x7ffUL
        );

    /*
     * Raw exponent 2047 means infinity or NaN.
     * Remap it to the largest finite exponent.
     */
    if (exponent == 0x7ff)
        exponent = 0x7fe;

    representation.exponent =
        exponent;

    representation.sign =
        (bits & (1UL << 63)) != 0;

    return representation.value;
}


P randomFinitePoint(ref ulong state)
    pure nothrow @safe @nogc
{
    return P(
        randomFiniteDouble(state),
        randomFiniteDouble(state),
        randomFiniteDouble(state)
    );
}


void main()
    @safe
{
    const P zero =
        P(
            0.0,
            0.0,
            0.0
        );


    /*
     * Canonical positive public convention:
     *
     * det(b-a, c-a, d-a) > 0.
     */
    checkExpected(
        +1,
        zero,
        P(1.0, 0.0, 0.0),
        P(0.0, 1.0, 0.0),
        P(0.0, 0.0, 1.0)
    );


    /*
     * Exact coplanarity.
     */
    checkExpected(
        0,
        zero,
        P(1.0, 0.0, 0.0),
        P(0.0, 1.0, 0.0),
        P(1.0, 1.0, 0.0)
    );


    /*
     * Complete degeneracy.
     */
    checkExpected(
        0,
        P(7.0, -2.0, 9.0),
        P(7.0, -2.0, 9.0),
        P(7.0, -2.0, 9.0),
        P(7.0, -2.0, 9.0)
    );


    /*
     * Smallest positive binary64 subnormal in all three basis
     * directions.
     *
     * Its exact determinant is positive although its magnitude is far
     * below what ordinary binary64 multiplication can represent.
     */
    enum double minSubnormal =
        0x0.0000000000001p-1022;

    checkExpected(
        +1,
        zero,
        P(minSubnormal, 0.0, 0.0),
        P(0.0, minSubnormal, 0.0),
        P(0.0, 0.0, minSubnormal)
    );


    /*
     * Maximum finite coordinates.
     *
     * The mathematical determinant magnitude is far beyond binary64.
     */
    checkExpected(
        +1,
        zero,
        P(double.max, 0.0, 0.0),
        P(0.0, double.max, 0.0),
        P(0.0, 0.0, double.max)
    );


    /*
     * Coordinate subtraction itself exceeds binary64:
     *
     *     double.max - (-double.max)
     *
     * Exact dyadic arithmetic remains valid.
     */
    checkExpected(
        +1,
        P(-double.max, 0.0, 0.0),
        P( double.max, 0.0, 0.0),
        P(-double.max, 1.0, 0.0),
        P(-double.max, 0.0, 1.0)
    );


    /*
     * Full-span exact coplanarity with overflowing X/Y differences.
     */
    checkExpected(
        0,
        P(-double.max, -double.max, 0.0),
        P( double.max, -double.max, 0.0),
        P(-double.max,  double.max, 0.0),
        P( double.max,  double.max, 0.0)
    );


    /*
     * Extreme exponent imbalance.
     */
    checkExpected(
        +1,
        zero,
        P(double.max, 0.0, 0.0),
        P(0.0, minSubnormal, 0.0),
        P(0.0, 0.0, double.max)
    );


    ulong state =
        0xa54f_f53a_5f1d_36f1UL;

    enum randomCases =
        5_000;


    /*
     * Arbitrary finite binary64 bit patterns.
     */
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
        5_000;


    /*
     * Arbitrary full-range X/Y coordinates constrained to z = 0.
     *
     * Every case is exactly coplanar irrespective of exponent range.
     */
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

        checkExpected(
            0,
            a,
            b,
            c,
            d
        );
    }


    writeln(
        "orientation3 full-range dyadic prototype: PASS (",
        randomCases,
        " arbitrary finite + ",
        coplanarCases,
        " full-range coplanar cases)"
    );
}
