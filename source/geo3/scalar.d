module geo3.scalar;

static import euclid_core.scalar;

alias isGeoScalar =
    euclid_core.scalar.isGeoScalar;

alias MetricScalar =
    euclid_core.scalar.MetricScalar;

alias IntersectionScalar =
    euclid_core.scalar.IntersectionScalar;

@safe unittest
{
    static assert(isGeoScalar!double);
    static assert(is(MetricScalar!double == double));
    static assert(is(IntersectionScalar!double == double));
}
