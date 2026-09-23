module app;

import geo3 : Point3;

import geo.internal.binary64_rounding :
    roundedSub;

import geo.internal.expansion :
    ExpansionBuffer,
    TwoComponent,
    expansionSign,
    fastExpansionSumZeroElim,
    negateExpansion,
    scaleExpansionZeroElim,
    twoDiff;

import std.bigint : BigInt;
import std.bitmanip : DoubleRep;
import std.math.traits : isFinite;
import std.stdio : writeln;


alias P = Point3!double;


/*
 * Conservative research working range.
 *
 * This stage is intentionally not the complete finite binary64 backend.
 * Inputs outside this range belong to the later exact dyadic fallback.
 */
enum double minWorkingMagnitude =
    0x1p-200;

enum double maxWorkingMagnitude =
    0x1p+200;


double magnitude(double value)
    pure nothrow @safe @nogc
{
    return value < 0.0
        ? -value
        : value;
}


bool supportedComponent(double value)
    pure nothrow @safe @nogc
{
    if (!isFinite(value))
        return false;

    if (value == 0.0)
        return true;

    const double absolute =
        magnitude(value);

    return absolute >= minWorkingMagnitude &&
           absolute <= maxWorkingMagnitude;
}


/*
 * Exact lhs-rhs expansion.
 *
 * Components are stored least-significant first.
 */
bool buildDifference(
    double lhs,
    double rhs,
    ref ExpansionBuffer!2 result
)
    pure nothrow @safe @nogc
{
    const double rounded =
        roundedSub(lhs, rhs);

    if (!isFinite(rounded))
        return false;

    const TwoComponent difference =
        twoDiff(lhs, rhs);

    if (
        !supportedComponent(difference.low) ||
        !supportedComponent(difference.high)
    )
        return false;

    result.clear();

    if (difference.low != 0.0)
        result.append(difference.low);

    if (difference.high != 0.0)
        result.append(difference.high);

    if (result.empty)
        result.append(0.0);

    return true;
}


/*
 * Exact multiplication by a two-component difference expansion.
 */
void multiplyByDifference(
    size_t SourceCapacity,
    size_t ResultCapacity
)(
    ref const ExpansionBuffer!SourceCapacity source,
    ref const ExpansionBuffer!2 difference,
    ref ExpansionBuffer!ResultCapacity result
)
    pure nothrow @safe @nogc
if (ResultCapacity >= 4 * SourceCapacity)
{
    assert(!source.empty);
    assert(!difference.empty);
    assert(difference.length <= 2);

    ExpansionBuffer!(2 * SourceCapacity) first;
    ExpansionBuffer!(2 * SourceCapacity) second;

    scaleExpansionZeroElim(
        source,
        difference[0],
        first
    );

    if (difference.length == 2)
    {
        scaleExpansionZeroElim(
            source,
            difference[1],
            second
        );
    }

    fastExpansionSumZeroElim(
        first,
        second,
        result
    );
}


/*
 * Exact product of three at-most-two-component coordinate differences.
 *
 * 2 x 2 -> capacity 8
 * 8 x 2 -> capacity 32
 */
void tripleProduct(
    ref const ExpansionBuffer!2 x,
    ref const ExpansionBuffer!2 y,
    ref const ExpansionBuffer!2 z,
    ref ExpansionBuffer!32 result
)
    pure nothrow @safe @nogc
{
    ExpansionBuffer!8 xy;

    multiplyByDifference(
        x,
        y,
        xy
    );

    multiplyByDifference(
        xy,
        z,
        result
    );
}


/*
 * Exact expansion evaluation of
 *
 *     det(b-a, c-a, d-a)
 *
 * within the conservative expansion working range.
 *
 * false means:
 *
 *     route this finite input to the later full-range dyadic backend.
 */
bool tryOrientationExactExpansion(
    P a,
    P b,
    P c,
    P d,
    out int sign
)
    pure nothrow @safe @nogc
{
    if (
        !a.isFinite ||
        !b.isFinite ||
        !c.isFinite ||
        !d.isFinite
    )
        return false;

    ExpansionBuffer!2 ux;
    ExpansionBuffer!2 uy;
    ExpansionBuffer!2 uz;

    ExpansionBuffer!2 vx;
    ExpansionBuffer!2 vy;
    ExpansionBuffer!2 vz;

    ExpansionBuffer!2 wx;
    ExpansionBuffer!2 wy;
    ExpansionBuffer!2 wz;

    if (
        !buildDifference(b.x, a.x, ux) ||
        !buildDifference(b.y, a.y, uy) ||
        !buildDifference(b.z, a.z, uz) ||
        !buildDifference(c.x, a.x, vx) ||
        !buildDifference(c.y, a.y, vy) ||
        !buildDifference(c.z, a.z, vz) ||
        !buildDifference(d.x, a.x, wx) ||
        !buildDifference(d.y, a.y, wy) ||
        !buildDifference(d.z, a.z, wz)
    )
        return false;


    ExpansionBuffer!32 t0;
    ExpansionBuffer!32 t1;
    ExpansionBuffer!32 t2;

    ExpansionBuffer!32 raw3;
    ExpansionBuffer!32 raw4;
    ExpansionBuffer!32 raw5;

    ExpansionBuffer!32 t3;
    ExpansionBuffer!32 t4;
    ExpansionBuffer!32 t5;


    tripleProduct(
        ux, vy, wz,
        t0
    );

    tripleProduct(
        uy, vz, wx,
        t1
    );

    tripleProduct(
        uz, vx, wy,
        t2
    );

    tripleProduct(
        uz, vy, wx,
        raw3
    );

    tripleProduct(
        uy, vx, wz,
        raw4
    );

    tripleProduct(
        ux, vz, wy,
        raw5
    );


    negateExpansion(
        raw3,
        t3
    );

    negateExpansion(
        raw4,
        t4
    );

    negateExpansion(
        raw5,
        t5
    );


    ExpansionBuffer!64 s01;
    ExpansionBuffer!64 s23;
    ExpansionBuffer!64 s45;

    fastExpansionSumZeroElim(
        t0,
        t1,
        s01
    );

    fastExpansionSumZeroElim(
        t2,
        t3,
        s23
    );

    fastExpansionSumZeroElim(
        t4,
        t5,
        s45
    );

    ExpansionBuffer!128 s0123;

    fastExpansionSumZeroElim(
        s01,
        s23,
        s0123
    );

    ExpansionBuffer!192 determinant;

    fastExpansionSumZeroElim(
        s0123,
        s45,
        determinant
    );

    sign =
        expansionSign(determinant);

    return true;
}


/*
 * Exact representation of one finite binary64 value as
 *
 *     integer * 2^-1074
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
 * Independent exact BigInt oracle for the geo3-d sign convention.
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


void requireExpansion(
    P a,
    P b,
    P c,
    P d
)
    @safe
{
    int sign;

    const bool success =
        tryOrientationExactExpansion(
            a,
            b,
            c,
            d,
            sign
        );

    assert(success);

    const int expected =
        oracleOrientation(
            a,
            b,
            c,
            d
        );

    assert(sign == expected);


    int reversed;

    assert(
        tryOrientationExactExpansion(
            a,
            c,
            b,
            d,
            reversed
        )
    );

    assert(reversed == -sign);
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


double randomModerate(ref ulong state)
    pure nothrow @safe @nogc
{
    return cast(double)
        cast(int)
        nextRandom(state);
}


P randomModeratePoint(ref ulong state)
    pure nothrow @safe @nogc
{
    return P(
        randomModerate(state),
        randomModerate(state),
        randomModerate(state)
    );
}


void main()
    @safe
{
    const P a =
        P(0.0, 0.0, 0.0);

    const P b =
        P(1.0, 0.0, 1.0);

    const P c =
        P(0.0, 1.0, 1.0);


    /*
     * Cancellation-heavy exact coplanarity.
     */
    requireExpansion(
        a,
        b,
        c,
        P(1.0, 1.0, 2.0)
    );


    /*
     * One ULP above the plane.
     */
    requireExpansion(
        a,
        b,
        c,
        P(
            1.0,
            1.0,
            0x1.0000000000001p+1
        )
    );


    /*
     * One representable value below 2.0.
     */
    requireExpansion(
        a,
        b,
        c,
        P(
            1.0,
            1.0,
            0x1.fffffffffffffp+0
        )
    );


    /*
     * Exercise a non-zero TwoDiff tail.
     */
    requireExpansion(
        P(0x1p-100, 0.0, 0.0),
        P(1.0, 0.0, 0.0),
        P(0.0, 1.0, 0.0),
        P(0.0, 0.0, 1.0)
    );


    /*
     * Complete degeneracy.
     */
    requireExpansion(
        P(7.0, -2.0, 9.0),
        P(7.0, -2.0, 9.0),
        P(7.0, -2.0, 9.0),
        P(7.0, -2.0, 9.0)
    );


    /*
     * Extreme finite values deliberately fall through.
     */
    int unusedSign;

    assert(
        !tryOrientationExactExpansion(
            P(-double.max, 0.0, 0.0),
            P( double.max, 0.0, 0.0),
            P(0.0, 1.0, 0.0),
            P(0.0, 0.0, 1.0),
            unusedSign
        )
    );


    enum double minSubnormal =
        0x0.0000000000001p-1022;

    assert(
        !tryOrientationExactExpansion(
            P(0.0, 0.0, 0.0),
            P(minSubnormal, 0.0, 0.0),
            P(0.0, minSubnormal, 0.0),
            P(0.0, 0.0, minSubnormal),
            unusedSign
        )
    );


    ulong state =
        0x3c6e_f372_fe94_f82bUL;

    enum randomCases =
        50_000;

    foreach (_; 0 .. randomCases)
    {
        requireExpansion(
            randomModeratePoint(state),
            randomModeratePoint(state),
            randomModeratePoint(state),
            randomModeratePoint(state)
        );
    }


    writeln(
        "orientation3 exact expansion prototype: PASS (",
        randomCases,
        " moderate random cases)"
    );
}
