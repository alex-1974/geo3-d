module app;

import geo3 :
    Orientation3,
    Point3,
    orientation;

import geo3.internal.orientation_dyadic :
    orientationDyadicExact;

import geo3.internal.orientation_expansion :
    tryOrientationExactExpansion;

import geo3.internal.orientation_filter :
    OrientationFilterResult,
    orientationFilter,
    orientationFilterFinite;

import geo3.internal.binary64_rounding :
    roundedAdd,
    roundedMul,
    roundedSub;

import bit_guard_filter :
    bitGuardFilter,
    bitGuardFilterFinite;

import candidate_orientation :
    candidateOrientation;

import production_filter_local;

import std.bitmanip :
    DoubleRep;

import core.time :
    MonoTime;

import std.stdio :
    writefln,
    writeln;


private struct Quad(T)
{
    Point3!T a;
    Point3!T b;
    Point3!T c;
    Point3!T d;
}


private ulong nextRandom(ref ulong state)
    pure nothrow @safe @nogc
{
    state ^= state << 13;
    state ^= state >> 7;
    state ^= state << 17;

    return state;
}


pragma(inline, false)
private int evalInt(
    ref const Quad!int q
)
    nothrow @safe @nogc
{
    return cast(int) orientation(
        q.a,
        q.b,
        q.c,
        q.d
    );
}


pragma(inline, false)
private int evalLong(
    ref const Quad!long q
)
    nothrow @safe @nogc
{
    return cast(int) orientation(
        q.a,
        q.b,
        q.c,
        q.d
    );
}


pragma(inline, false)
private int evalFloat(
    ref const Quad!float q
)
    nothrow @safe @nogc
{
    return cast(int) orientation(
        q.a,
        q.b,
        q.c,
        q.d
    );
}


pragma(inline, false)
private int evalDouble(
    ref const Quad!double q
)
    nothrow @safe @nogc
{
    return cast(int) orientation(
        q.a,
        q.b,
        q.c,
        q.d
    );
}


pragma(inline, false)
private int evalCandidateDouble(
    ref const Quad!double q
)
    nothrow @safe @nogc
{
    return cast(int) candidateOrientation(
        q.a,
        q.b,
        q.c,
        q.d
    );
}


pragma(inline, false)
private int evalFilter(
    ref const Quad!double q
)
    nothrow @safe @nogc
{
    return cast(int) orientationFilter(
        q.a,
        q.b,
        q.c,
        q.d
    );
}


pragma(inline, false)
private int evalLocalProductionFiniteFilter(
    ref const Quad!double q
)
    nothrow @safe @nogc
{
    return cast(int)
        production_filter_local.orientationFilterFinite(
            q.a,
            q.b,
            q.c,
            q.d
        );
}


pragma(inline, false)
private int evalProductionFiniteFilter(
    ref const Quad!double q
)
    nothrow @safe @nogc
{
    return cast(int) orientationFilterFinite(
        q.a,
        q.b,
        q.c,
        q.d
    );
}


pragma(inline, false)
private int evalBitGuardFilter(
    ref const Quad!double q
)
    nothrow @safe @nogc
{
    return cast(int) bitGuardFilter(
        q.a,
        q.b,
        q.c,
        q.d
    );
}


pragma(inline, false)
private int evalBitGuardFiniteFilter(
    ref const Quad!double q
)
    nothrow @safe @nogc
{
    return cast(int) bitGuardFilterFinite(
        q.a,
        q.b,
        q.c,
        q.d
    );
}


pragma(inline, false)
private int evalExpansion(
    ref const Quad!double q
)
    nothrow @safe @nogc
{
    int sign;

    const bool success =
        tryOrientationExactExpansion(
            q.a,
            q.b,
            q.c,
            q.d,
            sign
        );

    return success
        ? sign
        : 17;
}


pragma(inline, false)
private int evalDyadic(
    ref const Quad!double q
)
    nothrow @safe @nogc
{
    return orientationDyadicExact(
        q.a,
        q.b,
        q.c,
        q.d
    );
}


private void runBenchmark(
    alias evaluate,
    Q
)(
    string label,
    Q[] cases,
    size_t rounds
)
{
    long warmupChecksum;

    foreach (ref const q; cases)
    {
        warmupChecksum +=
            evaluate(q);
    }


    long checksum =
        warmupChecksum;

    const auto started =
        MonoTime.currTime;


    foreach (_; 0 .. rounds)
    {
        foreach (ref const q; cases)
        {
            checksum +=
                evaluate(q);
        }
    }


    const auto elapsed =
        MonoTime.currTime -
        started;

    const ulong calls =
        cast(ulong) cases.length *
        cast(ulong) rounds;

    const double nanoseconds =
        cast(double)
        elapsed.total!"nsecs";

    const double perCall =
        nanoseconds /
        cast(double) calls;


    writefln(
        "%-30s %12.3f ns/call  calls=%s  checksum=%s",
        label,
        perCall,
        calls,
        checksum
    );
}


private double randomFiniteDouble(
    ref ulong state
)
    pure nothrow @safe @nogc
{
    const ulong bits =
        nextRandom(state);


    DoubleRep representation;

    representation.value =
        0.0;


    representation.fraction =
        bits &
        (
            (1UL << 52) -
            1
        );


    ushort exponent =
        cast(ushort)(
            (
                bits >>
                52
            ) &
            0x7ffUL
        );


    if (exponent == 0x7ff)
    {
        exponent =
            0x7fe;
    }


    representation.exponent =
        exponent;

    representation.sign =
        (
            bits &
            (1UL << 63)
        ) != 0;


    return
        representation.value;
}


private Quad!double randomFiniteQuad(
    ref ulong state
)
    pure nothrow @safe @nogc
{
    alias P =
        Point3!double;


    return Quad!double(
        P(
            randomFiniteDouble(state),
            randomFiniteDouble(state),
            randomFiniteDouble(state)
        ),
        P(
            randomFiniteDouble(state),
            randomFiniteDouble(state),
            randomFiniteDouble(state)
        ),
        P(
            randomFiniteDouble(state),
            randomFiniteDouble(state),
            randomFiniteDouble(state)
        ),
        P(
            randomFiniteDouble(state),
            randomFiniteDouble(state),
            randomFiniteDouble(state)
        )
    );
}


private Quad!double nearCoplanarCase(
    bool reverse
)
{
    alias P =
        Point3!double;

    enum double aboveTwo =
        0x1.0000000000001p+1;


    auto result =
        Quad!double(
            P(
                0.0,
                0.0,
                0.0
            ),
            P(
                1.0,
                0.0,
                1.0
            ),
            P(
                0.0,
                1.0,
                1.0
            ),
            P(
                1.0,
                1.0,
                aboveTwo
            )
        );


    if (reverse)
    {
        const P temporary =
            result.b;

        result.b =
            result.c;

        result.c =
            temporary;
    }


    return result;
}


private Quad!double extremeCase(
    bool reverse
)
{
    alias P =
        Point3!double;


    auto result =
        Quad!double(
            P(
                -double.max,
                0.0,
                0.0
            ),
            P(
                0.0,
                1.0,
                0.0
            ),
            P(
                0.0,
                0.0,
                1.0
            ),
            P(
                double.max,
                0.0,
                0.0
            )
        );


    if (reverse)
    {
        const P temporary =
            result.b;

        result.b =
            result.c;

        result.c =
            temporary;
    }


    return result;
}



/*
 * EXPERIMENT ONLY.
 *
 * Same first-stage determinant and Shewchuk A-bound arithmetic as the
 * production filter, using the same explicit binary64 rounding operations.
 *
 * Deliberately omits the production overflow / underflow / subnormal guards.
 * Therefore this function is valid here only for the moderate benchmark
 * dataset. It must not be used as a production predicate.
 *
 * Its purpose is to measure the performance ceiling available from reducing
 * safety-check overhead while preserving the actual rounded arithmetic.
 */
private int uncheckedModerateFilter(
    ref const Quad!double q
)
    nothrow @safe @nogc
{
    enum double errboundA =
        0x1.c000000000007p-51;


    const double adx =
        roundedSub(q.a.x, q.d.x);

    const double bdx =
        roundedSub(q.b.x, q.d.x);

    const double cdx =
        roundedSub(q.c.x, q.d.x);

    const double ady =
        roundedSub(q.a.y, q.d.y);

    const double bdy =
        roundedSub(q.b.y, q.d.y);

    const double cdy =
        roundedSub(q.c.y, q.d.y);

    const double adz =
        roundedSub(q.a.z, q.d.z);

    const double bdz =
        roundedSub(q.b.z, q.d.z);

    const double cdz =
        roundedSub(q.c.z, q.d.z);


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


    const double det =
        roundedAdd(
            roundedAdd(
                termA,
                termB
            ),
            termC
        );


    const double absBdxcdy =
        bdxcdy < 0.0
            ? -bdxcdy
            : bdxcdy;

    const double absCdxbdy =
        cdxbdy < 0.0
            ? -cdxbdy
            : cdxbdy;

    const double absCdxady =
        cdxady < 0.0
            ? -cdxady
            : cdxady;

    const double absAdxcdy =
        adxcdy < 0.0
            ? -adxcdy
            : adxcdy;

    const double absAdxbdy =
        adxbdy < 0.0
            ? -adxbdy
            : adxbdy;

    const double absBdxady =
        bdxady < 0.0
            ? -bdxady
            : bdxady;

    const double absAdz =
        adz < 0.0
            ? -adz
            : adz;

    const double absBdz =
        bdz < 0.0
            ? -bdz
            : bdz;

    const double absCdz =
        cdz < 0.0
            ? -cdz
            : cdz;


    const double permanentA =
        roundedMul(
            roundedAdd(
                absBdxcdy,
                absCdxbdy
            ),
            absAdz
        );

    const double permanentB =
        roundedMul(
            roundedAdd(
                absCdxady,
                absAdxcdy
            ),
            absBdz
        );

    const double permanentC =
        roundedMul(
            roundedAdd(
                absAdxbdy,
                absBdxady
            ),
            absCdz
        );


    const double permanent =
        roundedAdd(
            roundedAdd(
                permanentA,
                permanentB
            ),
            permanentC
        );


    if (permanent == 0.0)
        return 0;


    const double errbound =
        roundedMul(
            errboundA,
            permanent
        );


    /*
     * Same sign normalization as production orientationFilter.
     */
    if (det > errbound)
        return -1;

    if (-det > errbound)
        return 1;

    return 2;
}


pragma(inline, false)
private int evalUncheckedModerateFilter(
    ref const Quad!double q
)
    nothrow @safe @nogc
{
    return uncheckedModerateFilter(q);
}


private void printPathDiagnostics()
{
    const auto near =
        nearCoplanarCase(false);

    const auto extreme =
        extremeCase(false);


    int nearSign;
    int extremeSign;


    const auto nearFilter =
        orientationFilter(
            near.a,
            near.b,
            near.c,
            near.d
        );

    const bool nearExpansion =
        tryOrientationExactExpansion(
            near.a,
            near.b,
            near.c,
            near.d,
            nearSign
        );


    const auto extremeFilter =
        orientationFilter(
            extreme.a,
            extreme.b,
            extreme.c,
            extreme.d
        );

    const bool extremeExpansion =
        tryOrientationExactExpansion(
            extreme.a,
            extreme.b,
            extreme.c,
            extreme.d,
            extremeSign
        );


    writeln(
        "=== PATH DIAGNOSTICS ==="
    );

    writefln(
        "near:    filter=%s expansion=%s expansion-sign=%s public=%s",
        cast(int) nearFilter,
        nearExpansion,
        nearSign,
        cast(int) orientation(
            near.a,
            near.b,
            near.c,
            near.d
        )
    );

    writefln(
        "extreme: filter=%s expansion=%s expansion-sign=%s public=%s dyadic=%s",
        cast(int) extremeFilter,
        extremeExpansion,
        extremeSign,
        cast(int) orientation(
            extreme.a,
            extreme.b,
            extreme.c,
            extreme.d
        ),
        orientationDyadicExact(
            extreme.a,
            extreme.b,
            extreme.c,
            extreme.d
        )
    );

    writeln();
}


void main()
{
    enum size_t regularCount =
        4096;

    enum size_t fallbackCount =
        256;


    auto intCases =
        new Quad!int[regularCount];

    auto longCases =
        new Quad!long[regularCount];

    auto floatCases =
        new Quad!float[regularCount];

    auto doubleCases =
        new Quad!double[regularCount];

    auto nearCases =
        new Quad!double[fallbackCount];

    auto extremeCases =
        new Quad!double[fallbackCount];


    ulong state =
        0x243f_6a88_85a3_08d3UL;


    foreach (index; 0 .. regularCount)
    {
        alias PI =
            Point3!int;

        alias PL =
            Point3!long;

        alias PF =
            Point3!float;

        alias PD =
            Point3!double;


        intCases[index] =
            Quad!int(
                PI(
                    cast(int) nextRandom(state),
                    cast(int) nextRandom(state),
                    cast(int) nextRandom(state)
                ),
                PI(
                    cast(int) nextRandom(state),
                    cast(int) nextRandom(state),
                    cast(int) nextRandom(state)
                ),
                PI(
                    cast(int) nextRandom(state),
                    cast(int) nextRandom(state),
                    cast(int) nextRandom(state)
                ),
                PI(
                    cast(int) nextRandom(state),
                    cast(int) nextRandom(state),
                    cast(int) nextRandom(state)
                )
            );


        longCases[index] =
            Quad!long(
                PL(
                    cast(long) nextRandom(state),
                    cast(long) nextRandom(state),
                    cast(long) nextRandom(state)
                ),
                PL(
                    cast(long) nextRandom(state),
                    cast(long) nextRandom(state),
                    cast(long) nextRandom(state)
                ),
                PL(
                    cast(long) nextRandom(state),
                    cast(long) nextRandom(state),
                    cast(long) nextRandom(state)
                ),
                PL(
                    cast(long) nextRandom(state),
                    cast(long) nextRandom(state),
                    cast(long) nextRandom(state)
                )
            );


        const int fax =
            cast(int) nextRandom(state);

        const int fay =
            cast(int) nextRandom(state);

        const int faz =
            cast(int) nextRandom(state);

        const int fbx =
            cast(int) nextRandom(state);

        const int fby =
            cast(int) nextRandom(state);

        const int fbz =
            cast(int) nextRandom(state);

        const int fcx =
            cast(int) nextRandom(state);

        const int fcy =
            cast(int) nextRandom(state);

        const int fcz =
            cast(int) nextRandom(state);

        const int fdx =
            cast(int) nextRandom(state);

        const int fdy =
            cast(int) nextRandom(state);

        const int fdz =
            cast(int) nextRandom(state);


        floatCases[index] =
            Quad!float(
                PF(
                    cast(float) fax,
                    cast(float) fay,
                    cast(float) faz
                ),
                PF(
                    cast(float) fbx,
                    cast(float) fby,
                    cast(float) fbz
                ),
                PF(
                    cast(float) fcx,
                    cast(float) fcy,
                    cast(float) fcz
                ),
                PF(
                    cast(float) fdx,
                    cast(float) fdy,
                    cast(float) fdz
                )
            );


        doubleCases[index] =
            Quad!double(
                PD(
                    cast(double) fax,
                    cast(double) fay,
                    cast(double) faz
                ),
                PD(
                    cast(double) fbx,
                    cast(double) fby,
                    cast(double) fbz
                ),
                PD(
                    cast(double) fcx,
                    cast(double) fcy,
                    cast(double) fcz
                ),
                PD(
                    cast(double) fdx,
                    cast(double) fdy,
                    cast(double) fdz
                )
            );
    }


    foreach (index; 0 .. fallbackCount)
    {
        nearCases[index] =
            nearCoplanarCase(
                (index & 1) != 0
            );

        extremeCases[index] =
            extremeCase(
                (index & 1) != 0
            );
    }


    printPathDiagnostics();


    /*
     * The unchecked prototype is benchmarked only on the moderate dataset.
     * Prove for every benchmark case that it produces exactly the same filter
     * classification as the production implementation.
     */
    foreach (ref const q; doubleCases)
    {
        assert(
            evalUncheckedModerateFilter(q) ==
            evalFilter(q)
        );
    }

    writeln(
        "moderate unchecked-filter equivalence: PASS"
    );


    foreach (ref const q; doubleCases)
    {
        assert(
            bitGuardFilter(
                q.a,
                q.b,
                q.c,
                q.d
            ) ==
            orientationFilter(
                q.a,
                q.b,
                q.c,
                q.d
            )
        );

        assert(
            bitGuardFilterFinite(
                q.a,
                q.b,
                q.c,
                q.d
            ) ==
            orientationFilter(
                q.a,
                q.b,
                q.c,
                q.d
            )
        );
    }


    foreach (ref const q; nearCases)
    {
        assert(
            bitGuardFilter(
                q.a,
                q.b,
                q.c,
                q.d
            ) ==
            orientationFilter(
                q.a,
                q.b,
                q.c,
                q.d
            )
        );

        assert(
            bitGuardFilterFinite(
                q.a,
                q.b,
                q.c,
                q.d
            ) ==
            orientationFilter(
                q.a,
                q.b,
                q.c,
                q.d
            )
        );
    }


    foreach (ref const q; extremeCases)
    {
        assert(
            bitGuardFilter(
                q.a,
                q.b,
                q.c,
                q.d
            ) ==
            orientationFilter(
                q.a,
                q.b,
                q.c,
                q.d
            )
        );

        assert(
            bitGuardFilterFinite(
                q.a,
                q.b,
                q.c,
                q.d
            ) ==
            orientationFilter(
                q.a,
                q.b,
                q.c,
                q.d
            )
        );
    }


    ulong validationState =
        0x1319_8a2e_0370_7344UL;

    enum size_t fullRangeValidationCases =
        100_000;


    foreach (_; 0 .. fullRangeValidationCases)
    {
        const auto q =
            randomFiniteQuad(
                validationState
            );


        const auto production =
            orientationFilter(
                q.a,
                q.b,
                q.c,
                q.d
            );


        assert(
            bitGuardFilter(
                q.a,
                q.b,
                q.c,
                q.d
            ) ==
            production
        );


        assert(
            bitGuardFilterFinite(
                q.a,
                q.b,
                q.c,
                q.d
            ) ==
            production
        );
    }


    alias PD =
        Point3!double;

    const auto nonFinite =
        Quad!double(
            PD(
                double.nan,
                0.0,
                0.0
            ),
            PD(
                1.0,
                0.0,
                0.0
            ),
            PD(
                0.0,
                1.0,
                0.0
            ),
            PD(
                0.0,
                0.0,
                1.0
            )
        );


    assert(
        bitGuardFilter(
            nonFinite.a,
            nonFinite.b,
            nonFinite.c,
            nonFinite.d
        ) ==
        orientationFilter(
            nonFinite.a,
            nonFinite.b,
            nonFinite.c,
            nonFinite.d
        )
    );


    writeln(
        "bit-guard equivalence: PASS"
    );

    writeln(
        "full-range finite cases: ",
        fullRangeValidationCases
    );


    /*
     * End-to-end public-equivalent validation.
     */
    foreach (ref const q; doubleCases)
    {
        assert(
            candidateOrientation(
                q.a,
                q.b,
                q.c,
                q.d
            ) ==
            orientation(
                q.a,
                q.b,
                q.c,
                q.d
            )
        );
    }


    foreach (ref const q; nearCases)
    {
        assert(
            candidateOrientation(
                q.a,
                q.b,
                q.c,
                q.d
            ) ==
            orientation(
                q.a,
                q.b,
                q.c,
                q.d
            )
        );
    }


    foreach (ref const q; extremeCases)
    {
        assert(
            candidateOrientation(
                q.a,
                q.b,
                q.c,
                q.d
            ) ==
            orientation(
                q.a,
                q.b,
                q.c,
                q.d
            )
        );
    }


    ulong pipelineValidationState =
        0xa409_3822_299f_31d0UL;

    enum size_t pipelineValidationCases =
        250_000;


    foreach (_; 0 .. pipelineValidationCases)
    {
        const auto q =
            randomFiniteQuad(
                pipelineValidationState
            );


        assert(
            candidateOrientation(
                q.a,
                q.b,
                q.c,
                q.d
            ) ==
            orientation(
                q.a,
                q.b,
                q.c,
                q.d
            )
        );
    }


    writeln(
        "candidate pipeline equivalence: PASS"
    );

    writeln(
        "candidate full-range cases: ",
        pipelineValidationCases
    );

    writeln();


    writeln(
        "=== PUBLIC API ==="
    );

    runBenchmark!evalInt(
        "int public",
        intCases,
        500
    );

    runBenchmark!evalLong(
        "long public",
        longCases,
        500
    );

    runBenchmark!evalFloat(
        "float public",
        floatCases,
        1_000
    );

    runBenchmark!evalDouble(
        "double filter public",
        doubleCases,
        1_000
    );

    runBenchmark!evalCandidateDouble(
        "double filter candidate",
        doubleCases,
        1_000
    );

    runBenchmark!evalDouble(
        "double near public",
        nearCases,
        1_000
    );

    runBenchmark!evalCandidateDouble(
        "double near candidate",
        nearCases,
        1_000
    );

    runBenchmark!evalDouble(
        "double extreme public",
        extremeCases,
        100
    );

    runBenchmark!evalCandidateDouble(
        "double extreme candidate",
        extremeCases,
        100
    );


    writeln();
    writeln(
        "=== INTERNAL STAGES ==="
    );

    runBenchmark!evalFilter(
        "filter direct",
        doubleCases,
        1_000
    );

    runBenchmark!evalProductionFiniteFilter(
        "filter production finite",
        doubleCases,
        1_000
    );

    runBenchmark!evalLocalProductionFiniteFilter(
        "filter production local",
        doubleCases,
        1_000
    );

    runBenchmark!evalUncheckedModerateFilter(
        "filter unchecked moderate",
        doubleCases,
        1_000
    );

    runBenchmark!evalBitGuardFilter(
        "filter bit guards",
        doubleCases,
        1_000
    );

    runBenchmark!evalBitGuardFiniteFilter(
        "filter bit guards finite",
        doubleCases,
        1_000
    );

    runBenchmark!evalExpansion(
        "expansion direct",
        nearCases,
        1_000
    );

    runBenchmark!evalDyadic(
        "dyadic direct",
        extremeCases,
        100
    );


    writeln();
    writeln(
        "=== FILTER TARGET-BOUNDARY TEST ==="
    );

    runBenchmark!evalProductionFiniteFilter(
        "TARGET library production 1",
        doubleCases,
        1_000
    );

    runBenchmark!evalLocalProductionFiniteFilter(
        "TARGET local production 1",
        doubleCases,
        1_000
    );

    runBenchmark!evalBitGuardFiniteFilter(
        "TARGET experiment 1",
        doubleCases,
        1_000
    );

    runBenchmark!evalBitGuardFiniteFilter(
        "TARGET experiment 2",
        doubleCases,
        1_000
    );

    runBenchmark!evalLocalProductionFiniteFilter(
        "TARGET local production 2",
        doubleCases,
        1_000
    );

    runBenchmark!evalProductionFiniteFilter(
        "TARGET library production 2",
        doubleCases,
        1_000
    );


    writeln();
    writeln(
        "=== DOUBLE FASTPATH ABBA ==="
    );

    runBenchmark!evalDouble(
        "ABBA A public 1",
        doubleCases,
        1_000
    );

    runBenchmark!evalCandidateDouble(
        "ABBA B candidate 1",
        doubleCases,
        1_000
    );

    runBenchmark!evalCandidateDouble(
        "ABBA B candidate 2",
        doubleCases,
        1_000
    );

    runBenchmark!evalDouble(
        "ABBA A public 2",
        doubleCases,
        1_000
    );


    writeln();
    writeln(
        "=== FILTER ABBA ==="
    );

    runBenchmark!evalProductionFiniteFilter(
        "ABBA A production finite 1",
        doubleCases,
        1_000
    );

    runBenchmark!evalBitGuardFiniteFilter(
        "ABBA B experiment finite 1",
        doubleCases,
        1_000
    );

    runBenchmark!evalBitGuardFiniteFilter(
        "ABBA B experiment finite 2",
        doubleCases,
        1_000
    );

    runBenchmark!evalProductionFiniteFilter(
        "ABBA A production finite 2",
        doubleCases,
        1_000
    );
}
