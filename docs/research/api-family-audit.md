# geo-d / geo3-d API family audit

## Status

Research document for the initial `geo3-d` API-family review.

This document does not authorize implementation merely for feature parity.
It classifies the relationship between the frozen `geo-d v2.0.0` API and the
current `geo3-d` skeleton.

The governing rule is:

> Symmetry where the mathematics is symmetric; specialization where it is not.

`geo-d` and `geo3-d` are independent sibling libraries, but API similarity
between equivalent concepts is a design contract.

## Reference baseline

2D reference:

- package: `geo-d`
- release: `v2.0.0`
- release commit:
  `83c974e0ca018bb378ce65a08ef9d3178decf9ad`
- public module root: `geo`

3D baseline at start of this audit:

- package: `geo3-d`
- branch: `main`
- commit:
  `aa64e3ae76106d3c7b15b905907a567414c31ad5`
- public module root: `geo3`

Shared contract package:

- `euclid-core-d ~>0.1.0`

## Classification vocabulary

The audit uses these classifications.

### aligned

The current 3D concept already follows the applicable family contract.

### natural 3D counterpart missing

The mathematics is dimensionally symmetric and the corresponding 3D API is
well defined, but the current `geo3-d` skeleton does not yet provide it.

This classification identifies family shape. It does not by itself authorize
implementation without the normal consumer/research gate.

### not aligned / repair required

The corresponding 3D API already exists, but its current behaviour does not
yet satisfy an applicable family contract and should be brought into
alignment unless subsequent research establishes a genuine mathematical
reason to diverge.

### intentional divergence

2D and 3D require deliberately different public semantics even though they
belong to one conceptual family.

### 2D-only

The current 2D concept does not imply a meaningful 3D counterpart.

### research required

A plausible 3D concept exists, but its mathematical, numerical, robustness,
or API semantics must be established before public API is introduced.

### consumer required

A plausible counterpart exists, but there is no current reason to add public
API until a concrete consumer requires it.

## Family invariants

Equivalent concepts should normally preserve the following across dimensions:

- scalar domain;
- exact value equality;
- `.init` semantics;
- finite/non-finite handling;
- explicit conversion policy;
- `MetricScalar` policy;
- `IntersectionScalar` policy where applicable;
- `try...` failure semantics;
- transactional output-parameter behaviour;
- allocation policy;
- ownership and view lifetime semantics;
- robustness guarantees;
- free-function naming;
- argument ordering;
- useful UFCS where natural;
- `pure`, `nothrow`, `@safe`, and `@nogc` where the equivalent 2D API provides
  them and 3D mathematics does not prevent them.

The dimension belongs in dimension-specific type names, not normally in
operation names.

## Structural audit

| Concept | geo-d v2.0.0 | geo3-d baseline | Status | Required family direction |
|---|---|---|---|---|
| point value representation | `Point2` | `Point3` | aligned | Preserve explicit `x/y/z`, zero `.init`, exact value equality, supported scalar domain |
| point finite-state query | `Point2.isFinite` | absent | natural 3D counterpart missing | `Point3.isFinite` should follow the same finite/non-finite policy over three coordinates |
| point affine algebra | `Point2` +/− `Vector2`, `Point2 - Point2 -> Vector2` | unavailable | natural 3D counterpart missing | Mirror the 2D affine contract once `Vector3` exists |
| forbidden point algebra | point + point, unary point negation, point scaling unavailable | currently unavailable | aligned | Preserve these exclusions when `Vector3` and affine point operations are added; add compile-negative regression tests |
| vector | `Vector2` | absent | natural 3D counterpart missing | `Vector3` should provide the corresponding vector-space algebra |
| segment | `Segment2` | absent | natural 3D counterpart missing | `Segment3` with endpoint order, valid degeneracy, exact equality and `isFinite` |
| axis-aligned bounds | `Bounds2` | absent | natural 3D counterpart missing | `Bounds3`; explicit empty state and component-wise XYZ bounds |
| polyline view | `Polyline2View` | absent | natural 3D counterpart missing | `Polyline3View`; same borrowing, indexing and segment traversal model |
| linear ring view | `LinearRing2View` | absent | natural 3D counterpart missing | `LinearRing3View` as a cyclic sequence of 3D points does not itself require planarity |
| polygon view | `Polygon2View` | absent | 2D-only | Do not infer `Polygon3View`; an embedded planar surface requires separate semantics and a consumer |
| polygon area | `signedArea`, `polygonArea`, `AreaScalar` | absent | 2D-only | No automatic 3D counterpart |
| point-in-polygon | `PointPolygonLocation`, `tryClassifyPointInPolygon` | absent | 2D-only | No automatic 3D counterpart |

## Metric audit

| Concept | geo-d v2.0.0 | geo3-d baseline | Status | Required family direction |
|---|---|---|---|---|
| `MetricScalar` | shared core identity | shared core identity | aligned | Keep one declaration identity through `euclid-core-d` |
| point `distance` | robust metric differencing policy + `hypot` | present | not aligned / repair required | API name/type shape aligns, but current numerical implementation does not yet match the 2D integer-difference policy |
| `squaredDistance` | present | absent | natural 3D counterpart missing | Same metric scalar and non-exact metric semantics |
| `segmentLength` | present | absent | natural 3D counterpart missing | Same operation name and segment-first UFCS shape |
| `polylineLength` | present | absent | natural 3D counterpart missing | Same stored-order summation semantics unless research shows a reason to differ |
| `tryNearestPoint` | segment / point operation | absent | natural 3D counterpart missing | Segment-first family form is mathematically natural in 3D |
| `tryPointSegmentDistance` | canonical segment-first v2 form | absent | natural 3D counterpart missing | Preserve `tryPointSegmentDistance(segment, point, result)` |

### Current `distance` numerical gap

The current 3D implementation computes each component difference as:

```text
cast(MetricScalar, coordinate A)
-
cast(MetricScalar, coordinate B)
```

For integral coordinates, this differs materially from the `geo-d v2.0.0`
policy, which obtains the signed coordinate difference without signed overflow
before conversion to metric floating-point arithmetic.

For sufficiently large `long` coordinates, two distinct nearby integers can
round to the same floating-point value when converted before subtraction.
The current 3D implementation can therefore lose a small non-zero coordinate
difference.

The nested `hypot` magnitude computation is useful, but does not repair
information already lost during coordinate differencing.

Before `geo3-d` metric API is treated as family-aligned, integer differencing
must receive the same numerical review as the 2D implementation.

## Conversion audit

| Concept | geo-d v2.0.0 | geo3-d baseline | Status | Required family direction |
|---|---|---|---|---|
| `tryConvert` point | present | absent | natural 3D counterpart missing | Same checked per-coordinate conversion rules extended to Z |
| `tryConvert` vector | present | absent | natural 3D counterpart missing | Same scalar conversion policy |
| `tryConvert` segment | present | absent | natural 3D counterpart missing | Same transactional endpoint conversion |
| explicit rounding helpers | `rounded`, `floored`, `ceiled`, `truncated` | absent | consumer required | Review exact 2D overload set before introducing matching 3D forms |

No implicit quantisation should be introduced merely to simplify conversion.

## Bounds audit

`Bounds3` is a natural axis-aligned 3D counterpart of `Bounds2`.

Expected family semantics include:

- `Bounds3.init` is empty;
- empty state is explicit rather than encoded with coordinate sentinels;
- a degenerate min == max bounds is non-empty;
- min <= max is component-wise over X, Y and Z;
- floating non-empty bounds may contain infinities but not NaN;
- empty bounds are finite vacuously;
- construction and accumulation use transactional `try...` behaviour;
- no allocation is required.

Natural `tryBounds` families exist for:

- `Segment3`;
- `Polyline3View`;
- `LinearRing3View`.

A polygon overload is not implied because `Polygon3View` is not implied.

Status: natural 3D counterpart missing.

## Orientation audit

### 2D

The established affine 2D operation is:

```d
orientation(a, b, c)
```

and returns `Orientation2`.

### 3D

The plausible affine 3D family operation is:

```d
orientation(a, b, c, d)
```

and would classify the sign of oriented tetrahedral volume.

A three-point 3D cross product is not this affine orientation operation.

The exact public semantics of `Orientation3`, including names for the
negative, zero/coplanar, and positive states, are not yet fixed by family
symmetry alone.

Robust implementation is also dimension-specific:

- integer determinant growth differs from 2D;
- floating-point filters and exact fallback need 3D analysis;
- near-coplanar and extreme-scale inputs need independent oracle tests.

Status: research required.

No public `Orientation3` should be introduced merely by mechanically extending
the 2D enum.

## Intersection audit

Shared declarations already exported by `geo3-d`:

- `IntersectionScalar`
- `SegmentIntersectionKind`

These have common declaration identity through `euclid-core-d`.

This does not mean that `geo3-d` already has segment-intersection capability.

3D segment relationships include:

- disjoint skew segments;
- coplanar disjoint segments;
- endpoint contact;
- interior point intersection;
- collinear contact;
- positive-length collinear overlap;
- degenerate segments.

The common high-level result kind may remain usable, but robust predicate,
construction, degeneracy and representability semantics require dedicated
3D review.

Status for algorithms: research required.

## Ring audit

`LinearRing3View` as a representation is mathematically well defined without
requiring planarity: it is a cyclic ordered sequence of 3D vertices with an
implicit closing segment.

Therefore representation and validation must be considered separately.

### Representation

Status: natural 3D counterpart missing.

Expected family properties:

- non-owning read-only view;
- caller-owned contiguous backing storage;
- empty `.init`;
- empty and degenerate representations permitted;
- implicit closure;
- no normalization;
- `length`;
- `empty`;
- `segmentCount`;
- indexing by value;
- `Segment3` traversal.

### Validation

The shared declarations:

- `RingValidationIssue`
- `RingValidationResult`

already have common declaration identity.

However, `validateRing` cannot be copied mechanically from the current 2D
implementation because its intersection predicates are dimension-specific.

A 3D validation design must define and robustly handle at least:

- too few vertices;
- non-finite coordinates;
- zero-length edges;
- point self-intersection;
- positive-length self-overlap;
- valid skew non-intersection;
- the expected shared endpoint of adjacent ring edges.

Planarity must not be silently added as a validity requirement merely because
2D polygon rings are planar by construction.

Status: research required.

## Simplification audit

`douglasPeuckerWorkspaceSize` already has shared declaration identity.

The actual `trySimplifyDouglasPeuckerInto` algorithm is absent from `geo3-d`.

Douglas-Peucker over a 3D polyline is mathematically natural when point-to-
segment distance is Euclidean 3D distance.

The family shape can therefore remain:

```d
trySimplifyDouglasPeuckerInto(
    Polyline3View,
    tolerance,
    destination,
    workspace,
    written
)
```

with caller-owned destination and workspace.

Implementation should not precede the required `Polyline3View`, `Segment3`
and point-to-segment metric foundations.

Status: consumer required.

The mathematical family counterpart is already clear; the consumer gate applies
to implementing and publishing the operation, not to its classification as a
natural 3D polyline algorithm.

## Shared-core audit

The current shared contract set is:

- `isGeoScalar`
- `MetricScalar`
- `IntersectionScalar`
- `SegmentIntersectionKind`
- `RingValidationIssue`
- `RingValidationResult`
- `douglasPeuckerWorkspaceSize`

Current `geo3-d` aliases preserve the common declaration identity.

Status: aligned.

No additional declaration should move into `euclid-core-d` merely because
2D and 3D implementations contain similar source code.

## Root-package capability versus contract exports

The current `geo3` root exports:

- `Point3`;
- `distance`;
- the seven shared contracts.

Several shared declarations are therefore visible before the corresponding
3D geometry capability exists.

Examples:

- `SegmentIntersectionKind` exists, but `Segment3` and segment-intersection
  algorithms do not;
- `RingValidationIssue` and `RingValidationResult` exist, but
  `LinearRing3View` and `validateRing` do not;
- `douglasPeuckerWorkspaceSize` exists, but `Polyline3View` and
  `trySimplifyDouglasPeuckerInto` do not.

This is currently an architecture/coexistence skeleton, not evidence that
those 3D capabilities have already been designed.

Future documentation must make that distinction explicit.

## Tooling and package-parity audit

The current `geo3-d` skeleton also differs from the released `geo-d v2`
package setup.

Items to review before a stable `geo3-d` release include:

- D compiler/frontend minimum;
- DIP1000 compiler flags;
- package homepage;
- production package description rather than experimental skeleton wording;
- `.gitignore` treatment of `dub.selections.json` — resolved on this audit
  branch by aligning with the `geo-d` library policy;
- README depth;
- CHANGELOG;
- ROADMAP;
- CONTRIBUTING;
- SECURITY;
- getting-started documentation;
- API conventions;
- ADR structure;
- API-surface audit;
- durable sibling-coexistence consumer;
- release/freeze documentation.

These are engineering-family concerns, not reasons to add geometry API.

## Initial sequencing

This audit suggests the following dependency order for future work.

### Foundation family

1. `Vector3`
2. complete `Point3` family semantics
3. `Segment3`
4. `Bounds3`
5. `Polyline3View`
6. `LinearRing3View`

### Dimension-neutral operation families

After their prerequisite types exist:

- `squaredDistance`
- corrected/family-aligned `distance`
- `segmentLength`
- `polylineLength`
- `tryNearestPoint`
- `tryPointSegmentDistance`
- `tryBounds`
- `tryConvert`

Each addition remains subject to a concrete consumer or a research-backed
need rather than API completion for its own sake.

### Research-gated 3D algorithms

Separately investigate:

- `Orientation3` / robust 3D orientation;
- robust 3D segment relationships and intersection;
- `validateRing` in 3D;
- near-degenerate and extreme-scale numerical behaviour.

### Explicit non-goals

Do not infer from the 2D API:

- `Polygon3View`;
- 3D polygon area;
- point-in-3D-polygon API;
- polyhedra;
- meshes;
- triangulation;
- BSP/CSG;
- generic public N-dimensional geometry;
- spatial indexes.

Those require their own consumers and design decisions.

## Audit conclusion

The current `geo3-d` repository successfully proves the package relationship
and shared declaration identity with `geo-d`, but it is still an architecture
skeleton rather than a family-complete 3D geometry foundation.

The next design work should preserve the established `geo-d v2.0.0` contract
where the mathematics is dimensionally symmetric and isolate genuine 3D
mathematics behind explicit research decisions.

The first natural foundation gap is `Vector3`, because completing `Vector3`
also permits `Point3` to implement the already-established affine family
contract.
