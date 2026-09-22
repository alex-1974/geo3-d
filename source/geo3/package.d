module geo3;

public import geo3.intersection :
    SegmentIntersectionKind;

public import geo3.metric :
    distance;

public import geo3.point :
    Point3;

public import geo3.scalar :
    IntersectionScalar,
    MetricScalar,
    isGeoScalar;

public import geo3.simplification :
    douglasPeuckerWorkspaceSize;

public import geo3.topology_validation :
    RingValidationIssue,
    RingValidationResult;

public import geo3.vector :
    Vector3;
