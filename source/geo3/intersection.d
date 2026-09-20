module geo3.intersection;

static import euclid_core.intersection;

alias SegmentIntersectionKind =
    euclid_core.intersection.SegmentIntersectionKind;

@safe unittest
{
    SegmentIntersectionKind kind =
        SegmentIntersectionKind.point;

    assert(
        kind ==
        SegmentIntersectionKind.point
    );
}
