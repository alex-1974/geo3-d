module candidate_orientation;

import bit_guard_filter :
    bitGuardFilterFinite;

import geo3 :
    Orientation3,
    Point3;

import geo3.internal.orientation_dyadic :
    orientationDyadicExact;

import geo3.internal.orientation_expansion :
    tryOrientationExactExpansion;

import geo3.internal.orientation_filter :
    OrientationFilterResult;


/*
 * EXPERIMENT ONLY.
 *
 * End-to-end candidate for Point3!double Orientation3.
 *
 * The only changed stage is the first-stage filter:
 *
 *     production:
 *         orientationFilter
 *
 *     candidate:
 *         bitGuardFilterFinite
 *
 * Expansion and full-range dyadic fallback are the unchanged production
 * implementations.
 */
private Orientation3 fromSign(int sign)
    pure nothrow @safe @nogc
{
    if (sign > 0)
        return Orientation3.positive;

    if (sign < 0)
        return Orientation3.negative;

    return Orientation3.coplanar;
}


Orientation3 candidateOrientation(
    Point3!double a,
    Point3!double b,
    Point3!double c,
    Point3!double d
)
    pure nothrow @safe @nogc
{
    /*
     * Mirror the public production precondition.
     *
     * These assertions disappear from release builds.
     */
    assert(a.isFinite);
    assert(b.isFinite);
    assert(c.isFinite);
    assert(d.isFinite);


    const OrientationFilterResult filtered =
        bitGuardFilterFinite(
            a,
            b,
            c,
            d
        );


    final switch (filtered)
    {
        case OrientationFilterResult.negative:
            return Orientation3.negative;

        case OrientationFilterResult.coplanar:
            return Orientation3.coplanar;

        case OrientationFilterResult.positive:
            return Orientation3.positive;

        case OrientationFilterResult.uncertain:
            break;
    }


    int sign;


    if (
        tryOrientationExactExpansion(
            a,
            b,
            c,
            d,
            sign
        )
    )
    {
        return fromSign(sign);
    }


    return fromSign(
        orientationDyadicExact(
            a,
            b,
            c,
            d
        )
    );
}
