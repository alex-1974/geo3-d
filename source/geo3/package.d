module geo3;

public import geo3.bounding_box :
    tryBounds;

public import geo3.bounds :
    Bounds3;

public import geo3.convert :
    tryConvert;

public import geo3.intersection :
    SegmentIntersectionKind;

public import geo3.linear_ring_view :
    LinearRing3View;

public import geo3.metric :
    distance,
    polylineLength,
    segmentLength,
    squaredDistance;

public import geo3.point :
    Point3;

public import geo3.polyline_view :
    Polyline3View;

public import geo3.scalar :
    IntersectionScalar,
    MetricScalar,
    isGeoScalar;

public import geo3.segment :
    Segment3;

public import geo3.simplification :
    douglasPeuckerWorkspaceSize;

public import geo3.topology_validation :
    RingValidationIssue,
    RingValidationResult;

public import geo3.vector :
    Vector3;
