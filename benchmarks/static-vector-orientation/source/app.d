module app;

import geo3.internal.orientation_expansion :
    tryOrientationExactExpansion;
import geo3.point :
    Point3;

import std.conv : to;
import std.stdio : stderr, writeln;

enum size_t caseCount = 256;

private struct Case
{
    Point3!double a;
    Point3!double b;
    Point3!double c;
    Point3!double d;
}

private void fillCases(ref Case[caseCount] cases)
    pure nothrow @safe @nogc
{
    foreach (i, ref item; cases)
    {
        const double t =
            cast(double)(i & 31);

        double dZ;

        final switch (i % 3)
        {
            case 0:
                dZ = 1.2;
                break;

            case 1:
                dZ = 1.2000000000001;
                break;

            case 2:
                dZ = -0.7;
                break;
        }

        item =
            Case(
                Point3!double(
                    t + 0.1,
                    t * 0.5 + 0.2,
                    -t * 0.25 + 0.3
                ),
                Point3!double(
                    t + 1.7,
                    t + 2.1,
                    -t * 0.25 + 0.9
                ),
                Point3!double(
                    t - 1.3,
                    t + 0.8,
                    -t * 0.25 + 2.4
                ),
                Point3!double(
                    t + 0.6,
                    t - 1.4,
                    -t * 0.25 + dZ
                )
            );
    }
}

pragma(inline, false)
extern(C) ulong bench_expansion(
    scope const Case[] cases,
    size_t rounds)
    @safe @nogc nothrow
{
    ulong checksum =
        0xCBF2_9CE4_8422_2325UL;

    foreach (round; 0 .. rounds)
    {
        foreach (i, item; cases)
        {
            int sign;

            const success =
                tryOrientationExactExpansion(
                    item.a,
                    item.b,
                    item.c,
                    item.d,
                    sign
                );

            assert(success);

            checksum ^=
                cast(ulong)(sign + 2) +
                (cast(ulong) i + 1) *
                    0x9E37_79B9_7F4A_7C15UL +
                (cast(ulong) round + 1) *
                    0xD6E8_FEB8_6659_FD93UL;

            checksum *=
                0x0000_0100_0000_01B3UL;
        }
    }

    return checksum;
}

void main(string[] args)
{
    if (args.length != 2)
    {
        stderr.writeln(
            "usage: geo3-static-vector-orientation-probe <rounds>");
        return;
    }

    const rounds =
        to!size_t(args[1]);

    Case[caseCount] cases;
    fillCases(cases);

    const checksum =
        bench_expansion(
            cases[],
            rounds
        );

    writeln(checksum);
}
