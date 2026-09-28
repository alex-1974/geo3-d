/**
 * Research-only branch-vs-branch performance probe for the containers-d
 * StaticVector composition experiment in geo3-d.
 *
 * The same source is compiled against geo3-d develop and the research branch.
 */
module static_vector_consumer_probe;

import geo3.internal.expansion :
    ExpansionBuffer,
    fastExpansionSumZeroElim,
    scaleExpansionZeroElim;
import geo3.internal.orientation_expansion :
    tryOrientationExactExpansion;
import geo3.point : Point3;

import std.conv : to;
import std.stdio : stderr, writeln;

private struct OrientationCase
{
    Point3!double a;
    Point3!double b;
    Point3!double c;
    Point3!double d;
}

private __gshared ExpansionBuffer!2[2] scaleInputs;
private __gshared double[2] scaleScalars;

private __gshared ExpansionBuffer!4[2] sumLeft;
private __gshared ExpansionBuffer!4[2] sumRight;

private __gshared OrientationCase[2] orientationCoplanarCases;
private __gshared OrientationCase[2] orientationNearCases;

private void prepareCases()
{
    foreach (index; 0 .. 2)
        scaleInputs[index].clear();

    scaleInputs[0].append(0x1p-104);
    scaleInputs[0].append(1.0);
    scaleScalars[0] = 1.0 + 0x1p-27;

    scaleInputs[1].append(-0x1p-103);
    scaleInputs[1].append(2.0);
    scaleScalars[1] = -0.5 + 0x1p-28;

    foreach (index; 0 .. 2)
    {
        sumLeft[index].clear();
        sumRight[index].clear();
    }

    sumLeft[0].append(0x1p-156);
    sumLeft[0].append(-0x1p-104);
    sumLeft[0].append(0x1p-52);
    sumLeft[0].append(1.0);

    sumRight[0].append(-0x1p-155);
    sumRight[0].append(0x1p-103);
    sumRight[0].append(-0x1p-51);
    sumRight[0].append(2.0);

    sumLeft[1].append(-0x1p-158);
    sumLeft[1].append(0x1p-106);
    sumLeft[1].append(-0x1p-54);
    sumLeft[1].append(-1.5);

    sumRight[1].append(0x1p-157);
    sumRight[1].append(-0x1p-105);
    sumRight[1].append(0x1p-53);
    sumRight[1].append(-2.5);

    alias P = Point3!double;

    orientationCoplanarCases[0] =
        OrientationCase(
            P(0.0, 0.0, 0.0),
            P(1.0, 0.0, 1.0),
            P(0.0, 1.0, 1.0),
            P(1.0, 1.0, 2.0)
        );

    orientationCoplanarCases[1] =
        OrientationCase(
            P(0.0, 0.0, 0.0),
            P(0.0, 1.0, 1.0),
            P(1.0, 0.0, 1.0),
            P(1.0, 1.0, 2.0)
        );

    enum double aboveTwo =
        0x1.0000000000001p+1;

    orientationNearCases[0] =
        OrientationCase(
            P(0.0, 0.0, 0.0),
            P(1.0, 0.0, 1.0),
            P(0.0, 1.0, 1.0),
            P(1.0, 1.0, aboveTwo)
        );

    orientationNearCases[1] =
        OrientationCase(
            P(0.0, 0.0, 0.0),
            P(0.0, 1.0, 1.0),
            P(1.0, 0.0, 1.0),
            P(1.0, 1.0, aboveTwo)
        );
}

pragma(inline, false)
extern(C) ulong bench_scale2(size_t rounds)
{
    ulong checksum;

    foreach (i; 0 .. rounds)
    {
        const index = i & 1;

        ExpansionBuffer!4 result;

        scaleExpansionZeroElim(
            scaleInputs[index],
            scaleScalars[index],
            result
        );

        checksum +=
            cast(ulong) result.length +
            7UL * cast(ulong)(
                result[result.length - 1] != 0.0
            );
    }

    return checksum;
}

pragma(inline, false)
extern(C) ulong bench_sum4x4(size_t rounds)
{
    ulong checksum;

    foreach (i; 0 .. rounds)
    {
        const index = i & 1;

        ExpansionBuffer!8 result;

        fastExpansionSumZeroElim(
            sumLeft[index],
            sumRight[index],
            result
        );

        checksum +=
            cast(ulong) result.length +
            7UL * cast(ulong)(
                result[result.length - 1] != 0.0
            );
    }

    return checksum;
}

private ulong runOrientation(
    scope const OrientationCase[] cases,
    size_t rounds)
{
    ulong checksum;

    foreach (i; 0 .. rounds)
    {
        const value = cases[i & 1];

        int sign;

        const success =
            tryOrientationExactExpansion(
                value.a,
                value.b,
                value.c,
                value.d,
                sign
            );

        checksum +=
            cast(ulong) success +
            cast(ulong)(sign + 1) * 3UL;
    }

    return checksum;
}

pragma(inline, false)
extern(C) ulong bench_orientation_coplanar(size_t rounds)
{
    return runOrientation(
        orientationCoplanarCases[],
        rounds
    );
}

pragma(inline, false)
extern(C) ulong bench_orientation_near(size_t rounds)
{
    return runOrientation(
        orientationNearCases[],
        rounds
    );
}

void main(string[] args)
{
    if (args.length != 3)
    {
        stderr.writeln(
            "usage: static-vector-consumer-probe <scale2|sum4x4|orientation-coplanar|orientation-near> <rounds>");
        return;
    }

    prepareCases();

    const rounds = to!size_t(args[2]);
    ulong checksum;

    final switch (args[1])
    {
        case "scale2":
            checksum = bench_scale2(rounds);
            break;

        case "sum4x4":
            checksum = bench_sum4x4(rounds);
            break;

        case "orientation-coplanar":
            checksum =
                bench_orientation_coplanar(rounds);
            break;

        case "orientation-near":
            checksum =
                bench_orientation_near(rounds);
            break;
    }

    writeln(args[1], " ", checksum);
}
