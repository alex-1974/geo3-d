module app;

import geo3 : Point3;

import std.bigint : BigInt;
import std.bitmanip : DoubleRep;
import std.math.traits :
    isFinite,
    isSubnormal;
import std.stdio : writeln;


version (LDC)
{
    import ldc.llvmasm : __ir_pure;
}
else
{
    import core.math : toPrec;
}


alias P = Point3!double;


/*
 * Public-sign filter result for the research prototype.
 *
 * positive / negative refer to:
 *
 *     sign(det(b-a, c-a, d-a))
 */
enum FilterResult : byte
{
    negative = -1,
    coplanar =  0,
    positive =  1,
    uncertain = 2,
}


/*
 * Shewchuk orient3d first-stage error coefficient:
 *
 *     epsilon = 2^-53
 *     o3derrboundA = (7 + 56 * epsilon) * epsilon
 *
 * Exact binary64 value:
 *
 *     0x1.c000000000007p-51
 */
enum double o3dErrboundA =
    0x1.c000000000007p-51;


double roundedAdd(double lhs, double rhs)
    pure nothrow @safe @nogc
{
    version (LDC)
    {
        return __ir_pure!(
            `%r = fadd double %0, %1
             ret double %r`,
            double
        )(lhs, rhs);
    }
    else
    {
        return toPrec!double(lhs + rhs);
    }
}


double roundedSub(double lhs, double rhs)
    pure nothrow @safe @nogc
{
    version (LDC)
    {
        return __ir_pure!(
            `%r = fsub double %0, %1
             ret double %r`,
            double
        )(lhs, rhs);
    }
    else
    {
        return toPrec!double(lhs - rhs);
    }
}


double roundedMul(double lhs, double rhs)
    pure nothrow @safe @nogc
{
    version (LDC)
    {
        return __ir_pure!(
            `%r = fmul double %0, %1
             ret double %r`,
            double
        )(lhs, rhs);
    }
    else
    {
        return toPrec!double(lhs * rhs);
    }
}


double absolute(double value)
    pure nothrow @safe @nogc
{
    return value < 0.0
        ? -value
        : value;
}


bool normalOrZero(double value)
    pure nothrow @safe @nogc
{
    return value == 0.0 || (
        isFinite(value) &&
        !isSubnormal(value)
    );
}


bool productUnderflowed(
    double lhs,
    double rhs,
    double product
)
    pure nothrow @safe @nogc
{
    return product == 0.0 &&
           lhs != 0.0 &&
           rhs != 0.0;
}


FilterResult orientationFilter(
    P a,
    P b,
    P c,
    P d
)
    pure nothrow @safe @nogc
{
    if (
        !a.isFinite ||
        !b.isFinite ||
        !c.isFinite ||
        !d.isFinite
    )
        return FilterResult.uncertain;


    /*
     * Use the classical orient3d difference layout relative to d.
     *
     * Its determinant has the opposite sign from geo3-d's chosen
     * public convention:
     *
     *     det(b-a, c-a, d-a)
     *
     * Therefore certified signs are inverted before return.
     */
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
        !normalOrZero(adx) ||
        !normalOrZero(bdx) ||
        !normalOrZero(cdx) ||
        !normalOrZero(ady) ||
        !normalOrZero(bdy) ||
        !normalOrZero(cdy) ||
        !normalOrZero(adz) ||
        !normalOrZero(bdz) ||
        !normalOrZero(cdz)
    )
        return FilterResult.uncertain;


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
        productUnderflowed(bdx, cdy, bdxcdy) ||
        productUnderflowed(cdx, bdy, cdxbdy) ||
        productUnderflowed(cdx, ady, cdxady) ||
        productUnderflowed(adx, cdy, adxcdy) ||
        productUnderflowed(adx, bdy, adxbdy) ||
        productUnderflowed(bdx, ady, bdxady)
    )
        return FilterResult.uncertain;


    if (
        !normalOrZero(bdxcdy) ||
        !normalOrZero(cdxbdy) ||
        !normalOrZero(cdxady) ||
        !normalOrZero(adxcdy) ||
        !normalOrZero(adxbdy) ||
        !normalOrZero(bdxady)
    )
        return FilterResult.uncertain;


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
        !normalOrZero(bc) ||
        !normalOrZero(ca) ||
        !normalOrZero(ab)
    )
        return FilterResult.uncertain;


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
        productUnderflowed(adz, bc, termA) ||
        productUnderflowed(bdz, ca, termB) ||
        productUnderflowed(cdz, ab, termC)
    )
        return FilterResult.uncertain;


    if (
        !normalOrZero(termA) ||
        !normalOrZero(termB) ||
        !normalOrZero(termC)
    )
        return FilterResult.uncertain;


    const double detAB =
        roundedAdd(
            termA,
            termB
        );

    if (!normalOrZero(detAB))
        return FilterResult.uncertain;


    const double det =
        roundedAdd(
            detAB,
            termC
        );

    if (!normalOrZero(det))
        return FilterResult.uncertain;


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
        !normalOrZero(permanentPairA) ||
        !normalOrZero(permanentPairB) ||
        !normalOrZero(permanentPairC)
    )
        return FilterResult.uncertain;


    const double permanentA =
        roundedMul(
            permanentPairA,
            absolute(adz)
        );

    const double permanentB =
        roundedMul(
            permanentPairB,
            absolute(bdz)
        );

    const double permanentC =
        roundedMul(
            permanentPairC,
            absolute(cdz)
        );


    if (
        productUnderflowed(
            permanentPairA,
            absolute(adz),
            permanentA
        ) ||
        productUnderflowed(
            permanentPairB,
            absolute(bdz),
            permanentB
        ) ||
        productUnderflowed(
            permanentPairC,
            absolute(cdz),
            permanentC
        )
    )
        return FilterResult.uncertain;


    if (
        !normalOrZero(permanentA) ||
        !normalOrZero(permanentB) ||
        !normalOrZero(permanentC)
    )
        return FilterResult.uncertain;


    const double permanentAB =
        roundedAdd(
            permanentA,
            permanentB
        );

    if (!normalOrZero(permanentAB))
        return FilterResult.uncertain;


    const double permanent =
        roundedAdd(
            permanentAB,
            permanentC
        );

    if (!normalOrZero(permanent))
        return FilterResult.uncertain;


    /*
     * With every potentially underflowed multiplication rejected,
     * permanent == 0 implies that every exact determinant term is
     * structurally zero.
     */
    if (permanent == 0.0)
        return FilterResult.coplanar;


    const double errbound =
        roundedMul(
            o3dErrboundA,
            permanent
        );


    /*
     * Do not certify through overflow or an underflowed error bound.
     */
    if (
        !isFinite(errbound) ||
        errbound == 0.0 ||
        isSubnormal(errbound)
    )
        return FilterResult.uncertain;


    /*
     * det is the classical Shewchuk sign convention and is opposite
     * to geo3-d's chosen public determinant convention.
     */
    if (det > errbound)
        return FilterResult.negative;

    if (-det > errbound)
        return FilterResult.positive;

    return FilterResult.uncertain;
}


/*
 * Decode one finite binary64 value exactly as an integer multiple
 * of 2^-1074.
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
 * Independent exact oracle for the public convention:
 *
 *     sign(det(b-a, c-a, d-a))
 */
int oracleOrientation(
    P a,
    P b,
    P c,
    P d
)
    @safe
{
    assert(a.isFinite);
    assert(b.isFinite);
    assert(c.isFinite);
    assert(d.isFinite);

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


    const BigInt ux =
        bx - ax;

    const BigInt uy =
        by - ay;

    const BigInt uz =
        bz - az;

    const BigInt vx =
        cx - ax;

    const BigInt vy =
        cy - ay;

    const BigInt vz =
        cz - az;

    const BigInt wx =
        dx - ax;

    const BigInt wy =
        dy - ay;

    const BigInt wz =
        dz - az;


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


struct FilterCounts
{
    size_t negative;
    size_t coplanar;
    size_t positive;
    size_t uncertain;
}


void checkFilter(
    P a,
    P b,
    P c,
    P d,
    ref FilterCounts counts
)
    @safe
{
    const int exact =
        oracleOrientation(
            a,
            b,
            c,
            d
        );

    const FilterResult filtered =
        orientationFilter(
            a,
            b,
            c,
            d
        );


    final switch (filtered)
    {
        case FilterResult.negative:
            ++counts.negative;
            break;

        case FilterResult.coplanar:
            ++counts.coplanar;
            break;

        case FilterResult.positive:
            ++counts.positive;
            break;

        case FilterResult.uncertain:
            ++counts.uncertain;
            return;
    }


    /*
     * Central research invariant:
     *
     * Every certified filter result must equal the exact oracle.
     */
    assert(
        cast(int) filtered ==
        exact
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

    if (exponent == 0x7ff)
        exponent = 0x7fe;

    representation.exponent =
        exponent;

    representation.sign =
        (bits & (1UL << 63)) != 0;

    return representation.value;
}


double randomModerateDouble(ref ulong state)
    pure nothrow @safe @nogc
{
    return cast(double)
        cast(int)
        nextRandom(state);
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


P randomModeratePoint(ref ulong state)
    pure nothrow @safe @nogc
{
    return P(
        randomModerateDouble(state),
        randomModerateDouble(state),
        randomModerateDouble(state)
    );
}


void main()
    @safe
{
    /*
     * Canonical geo3-d positive convention.
     */
    const P a =
        P(0.0, 0.0, 0.0);

    const P b =
        P(1.0, 0.0, 0.0);

    const P c =
        P(0.0, 1.0, 0.0);

    const P d =
        P(0.0, 0.0, 1.0);

    assert(
        oracleOrientation(a, b, c, d) ==
        1
    );

    assert(
        orientationFilter(a, b, c, d) ==
        FilterResult.positive
    );

    assert(
        orientationFilter(a, c, b, d) ==
        FilterResult.negative
    );


    /*
     * Structurally exact coplanarity.
     */
    const P coplanarD =
        P(1.0, 1.0, 0.0);

    assert(
        oracleOrientation(
            a,
            b,
            c,
            coplanarD
        ) == 0
    );

    assert(
        orientationFilter(
            a,
            b,
            c,
            coplanarD
        ) == FilterResult.coplanar
    );


    /*
     * Smallest positive binary64 perturbation away from the plane.
     * The filter deliberately refuses subnormal arithmetic.
     */
    enum double minSubnormal =
        0x0.0000000000001p-1022;

    const P nearD =
        P(
            0.25,
            0.25,
            minSubnormal
        );

    assert(
        oracleOrientation(
            a,
            b,
            c,
            nearD
        ) == 1
    );

    assert(
        orientationFilter(
            a,
            b,
            c,
            nearD
        ) == FilterResult.uncertain
    );


    /*
     * Product overflow remains uncertain.
     */
    const P hugeB =
        P(
            double.max,
            0.0,
            0.0
        );

    const P hugeC =
        P(
            0.0,
            double.max,
            0.0
        );

    const P hugeD =
        P(
            0.0,
            0.0,
            double.max
        );

    assert(
        oracleOrientation(
            a,
            hugeB,
            hugeC,
            hugeD
        ) == 1
    );

    assert(
        orientationFilter(
            a,
            hugeB,
            hugeC,
            hugeD
        ) == FilterResult.uncertain
    );


    /*
     * Even coordinate subtraction may overflow although every input
     * coordinate is finite.
     */
    const P overflowA =
        P(
            -double.max,
            0.0,
            0.0
        );

    const P overflowB =
        P(
            0.0,
            1.0,
            0.0
        );

    const P overflowC =
        P(
            0.0,
            0.0,
            1.0
        );

    const P overflowD =
        P(
            double.max,
            0.0,
            0.0
        );

    assert(
        oracleOrientation(
            overflowA,
            overflowB,
            overflowC,
            overflowD
        ) == 1
    );

    assert(
        orientationFilter(
            overflowA,
            overflowB,
            overflowC,
            overflowD
        ) == FilterResult.uncertain
    );


    /*
     * Non-finite values are never certified.
     */
    assert(
        orientationFilter(
            P(double.nan, 0.0, 0.0),
            b,
            c,
            d
        ) == FilterResult.uncertain
    );

    assert(
        orientationFilter(
            P(double.infinity, 0.0, 0.0),
            b,
            c,
            d
        ) == FilterResult.uncertain
    );


    ulong state =
        0xbb67_ae85_84ca_a73bUL;

    FilterCounts moderateCounts;
    FilterCounts fullRangeCounts;
    FilterCounts coplanarCounts;

    enum moderateCases = 50_000;
    enum fullRangeCases = 5_000;
    enum coplanarCases = 10_000;


    /*
     * Ordinary-range inputs should normally be certified by the
     * first-stage filter.
     */
    foreach (_; 0 .. moderateCases)
    {
        checkFilter(
            randomModeratePoint(state),
            randomModeratePoint(state),
            randomModeratePoint(state),
            randomModeratePoint(state),
            moderateCounts
        );
    }


    /*
     * Arbitrary finite binary64 bit patterns exercise overflow,
     * underflow, exponent imbalance and the conservative uncertain path.
     */
    foreach (_; 0 .. fullRangeCases)
    {
        checkFilter(
            randomFinitePoint(state),
            randomFinitePoint(state),
            randomFinitePoint(state),
            randomFinitePoint(state),
            fullRangeCounts
        );
    }


    /*
     * Random exact z=0 geometry must always be exactly coplanar.
     */
    foreach (_; 0 .. coplanarCases)
    {
        const P ca =
            P(
                randomModerateDouble(state),
                randomModerateDouble(state),
                0.0
            );

        const P cb =
            P(
                randomModerateDouble(state),
                randomModerateDouble(state),
                0.0
            );

        const P cc =
            P(
                randomModerateDouble(state),
                randomModerateDouble(state),
                0.0
            );

        const P cd =
            P(
                randomModerateDouble(state),
                randomModerateDouble(state),
                0.0
            );

        checkFilter(
            ca,
            cb,
            cc,
            cd,
            coplanarCounts
        );
    }


    writeln(
        "orientation3 double filter prototype: PASS"
    );

    writeln(
        "moderate: positive=",
        moderateCounts.positive,
        " negative=",
        moderateCounts.negative,
        " coplanar=",
        moderateCounts.coplanar,
        " uncertain=",
        moderateCounts.uncertain
    );

    writeln(
        "full-range: positive=",
        fullRangeCounts.positive,
        " negative=",
        fullRangeCounts.negative,
        " coplanar=",
        fullRangeCounts.coplanar,
        " uncertain=",
        fullRangeCounts.uncertain
    );

    writeln(
        "coplanar: positive=",
        coplanarCounts.positive,
        " negative=",
        coplanarCounts.negative,
        " coplanar=",
        coplanarCounts.coplanar,
        " uncertain=",
        coplanarCounts.uncertain
    );
}
