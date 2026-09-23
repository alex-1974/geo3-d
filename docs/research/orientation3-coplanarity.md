# Orientation3 and coplanarity research

Status: research
Branch: research/orientation3-coplanarity
Baseline: geo3-d main 761dbc8e7bedb26ee1daa1a20e667cafae4d2cd3

## Purpose

This document investigates the public semantics and numerical architecture
required for a robust three-dimensional orientation predicate in geo3-d.

No public API described here is considered committed until the research
questions in this document are resolved.

The work is intentionally separated from:

- three-dimensional segment intersection;
- ring validation;
- planar polygon or surface abstractions;
- general plane types;
- mesh, tetrahedron, or polyhedron APIs.

Those consumers may later depend on Orientation3, but they must not dictate
premature additional public geometry types.

## Family context

geo-d v2.0.0 defines the robust two-dimensional affine predicate

    orientation(a, b, c)

returning Orientation2.

Its contract is a discrete predicate result, not an approximate determinant
magnitude.

The established numerical policy is:

- exact results for integral coordinates;
- certified robust results for supported floating-point coordinates;
- no global epsilon;
- no conversion of robust predicates into metric approximations;
- no allocation;
- pure, nothrow, @safe and @nogc where practical;
- real is not accepted by a robust predicate until a backend exists that
  preserves the required guarantee for the actual D real representation.

The corresponding three-dimensional affine predicate naturally requires four
points rather than three.

## Mathematical definition

For four points a, b, c and d, define

    D = det(b - a, c - a, d - a)

Equivalently,

    D = ((b - a) x (c - a)) . (d - a)

The sign of D is the three-dimensional affine orientation.

The magnitude of D is six times the signed tetrahedral volume, but the public
orientation predicate does not expose that magnitude.

The proposed semantic classification is:

    D < 0  -> negative
    D = 0  -> coplanar
    D > 0  -> positive

This sign convention must be part of the public contract.

## Relationship to established orient3d implementations

Shewchuk's classical ORIENT3D formulation evaluates the equivalent affine
determinant using a different ordering convention based on differences from
d.

For the canonical example

    a = (0, 0, 0)
    b = (1, 0, 0)
    c = (0, 1, 0)
    d = (0, 0, 1)

geo3-d's proposed determinant

    det(b - a, c - a, d - a)

is positive.

The classical Shewchuk determinant convention has the opposite sign for the
same ordered points.

Therefore an implementation adapted from Shewchuk must explicitly normalize
its sign to the geo3-d public convention.

The backend convention must never leak into the public API.

## Proposed public shape

The natural point API is

    Orientation3 orientation(
        Point3!T a,
        Point3!T b,
        Point3!T c,
        Point3!T d
    )

for supported robust scalar types.

No three-point Point3 overload is proposed.

Three points in 3D define, at most, an oriented normal direction or an
embedded plane orientation. They do not define the same affine orientation
predicate as four points.

A public Vector3 overload is also not proposed at this stage. Although the
scalar triple product gives a natural orientation of three vectors, geo-d's
existing family contract is point based and there is currently no concrete
consumer requiring a separate vector overload.

## Orientation3 result type

The preferred semantic names are:

- negative
- coplanar
- positive

with numerical values matching determinant sign:

    negative = -1
    coplanar = 0
    positive = 1

These names are preferable to left/right or above/below:

- left/right is intrinsically two-dimensional;
- above/below depends on a geometric viewing convention;
- positive/negative maps directly to the documented determinant.

The enum declaration should be ordered as:

    coplanar = 0
    negative = -1
    positive = 1

D defines the default initializer of a named enum as its first member.
Therefore:

    Orientation3.init == Orientation3.coplanar

This is deliberate.

Orientation2 preserves a historical non-neutral default for compatibility.
Orientation3 is a new type and has no corresponding compatibility constraint.
For a new determinant-sign classification, the zero/coplanar state is the
least surprising default.

Orientation3 is dimension-specific and should remain owned by geo3-d rather
than euclid-core-d.

## Coplanarity semantics

For the affine determinant,

    orientation(a, b, c, d) == Orientation3.coplanar

means that the four points are affinely dependent in 3D and therefore lie in
an affine subspace of dimension at most two.

This remains well defined for degenerate inputs.

Examples include:

- four ordinary coplanar points;
- three collinear points plus any fourth point;
- repeated points;
- all four points collinear;
- all four points equal.

The predicate therefore does not require a, b and c to define a unique plane.

This is different from a predicate that asks which side of a known,
non-degenerate oriented plane contains d.

No separate public coplanar(a, b, c, d) helper should be introduced in the
initial API.

The canonical spelling is:

    orientation(a, b, c, d) == Orientation3.coplanar

A later concrete consumer may justify a convenience predicate independently.

This must not be confused with an in-plane coplanar-orientation predicate.
Such a predicate answers a different question: orientation of points within
an already established plane. It is not an alias for the zero case of
Orientation3.

## Permutation identities

Orientation3 is alternating in its four point arguments.

Swapping any two arguments reverses the sign.

Every even permutation preserves the sign.

Every odd permutation reverses the sign.

This differs in an important way from the common three-argument Orientation2
cyclic identity.

For four arguments the cyclic rotation

    (a, b, c, d) -> (b, c, d, a)

is a four-cycle and therefore reverses the sign.

These identities should form part of the property-test suite.

## Translation invariance

Adding the same vector to all four points must preserve orientation exactly.

The implementation should preferably evaluate coordinate differences before
forming products so that this affine structure is explicit.

## Integral exactness

The complete int and long coordinate domains must produce mathematically
exact orientation classifications.

No signed source-coordinate subtraction may overflow.

### int

For signed 32-bit coordinates, one coordinate difference can require the full
32-bit unsigned magnitude.

A triple product of differences therefore requires fewer than 96 magnitude
bits.

Even the conservative sum of all six determinant-product magnitudes requires
fewer than 99 bits.

A private 128-bit exact accumulator is therefore sufficient for the complete
Point3!int domain.

### long

For signed 64-bit coordinates, one difference can require the full 64-bit
unsigned magnitude.

A triple product can require nearly 192 magnitude bits.

A conservative accumulation of all six determinant-product magnitudes
requires fewer than 195 bits.

A private 256-bit exact accumulator is therefore sufficient for the complete
Point3!long domain.

These are sufficient implementation bounds, not claims that the exact maximum
geometric determinant requires every one of those bits.

The implementation should prefer fixed-width exact arithmetic over production
BigInt allocation.

std.bigint.BigInt remains appropriate as an independent unittest oracle.

## Exact determinant backend design

The integral implementation should evaluate the determinant sign without
materializing a signed wide integer.

Let:

    u = b - a
    v = c - a
    w = d - a

and expand:

    det(u, v, w)
      = ux * vy * wz
      + uy * vz * wx
      + uz * vx * wy
      - uz * vy * wx
      - uy * vx * wz
      - ux * vz * wy

Each coordinate difference is represented as:

- a sign;
- an exact unsigned magnitude.

Each triple product therefore also has:

- a sign obtained from its three factor signs;
- an exact unsigned fixed-width magnitude.

Rather than summing signed wide integers, terms are accumulated into two
unsigned buckets:

    positiveSum
    negativeSum

The determinant sign is then:

- positive if positiveSum > negativeSum;
- negative if positiveSum < negativeSum;
- zero if they are equal.

This avoids requiring a signed arbitrary-width integer representation.

### Fixed-width limb model

A portable private fixed-width integer using 32-bit limbs is sufficient.

For non-zero component factors, the product of the signs of all six expanded
Leibniz terms is always negative:

- the product of the six permutation coefficients is -1;
- every component sign occurs in exactly two terms and therefore cancels
  from the total sign product.

Consequently the six non-zero terms cannot all have the same sign. Their
positive/negative split is necessarily 5/1, 3/3, or 1/5. Zero factors only
reduce the number of accumulated terms.

For Point3!int:

- one coordinate-difference magnitude is less than 2^32 and fits in 32 bits;
- one triple product is less than 2^96 and fits in 96 bits;
- either sign bucket receives at most five non-zero Leibniz terms;
- either bucket is therefore less than 5 * 2^96, which is less than 2^99;
- 99 magnitude bits are sufficient for a complete bucket;
- four 32-bit limbs provide 128 bits and are sufficient.

For Point3!long:

- one coordinate-difference magnitude is less than 2^64 and fits in 64 bits;
- one triple product is less than 2^192 and fits in 192 bits;
- either sign bucket receives at most five non-zero Leibniz terms;
- either bucket is therefore less than 5 * 2^192, which is less than 2^195;
- 195 magnitude bits are sufficient for a complete bucket;
- seven 32-bit limbs provide 224 bits and are sufficient.

Using eight limbs as a 256-bit implementation type would also be correct, but
is not mathematically required.

The existing geo-d internal UIntFixed implementation is a useful engineering
reference. It should not automatically be moved into euclid-core-d: this is
implementation machinery, not a declaration requiring common public identity.

geo3-d may initially provide its own private fixed-width implementation.
Implementation sharing can be reconsidered separately if repeated maintenance
becomes a concrete problem.

### Shared exact-sign engine

The determinant-sign machinery should be designed generically enough that the
same internal sign-plus-magnitude accumulation model can support:

- exact int orientation;
- exact long orientation;
- the full-range binary64 dyadic fallback.

This does not imply any public generic N-dimensional or arbitrary-precision
API.

### Full-range binary64 dyadic fallback

Every finite IEEE binary64 value is an exact integer multiple of 2^-1074.

When expressed in those common units:

- one finite coordinate requires at most 2098 magnitude bits;
- a difference between two finite coordinates requires at most 2099 bits;
- 66 32-bit limbs are therefore sufficient for an exact difference magnitude.

A triple product of three such differences is less than 2^6297 and may
therefore require 6297 magnitude bits.

The same Leibniz-sign invariant applies here: either sign bucket receives at
most five non-zero terms. A complete bucket is therefore less than

    5 * 2^6297

which is less than 2^6300.

Thus 6300 magnitude bits are sufficient for a complete exact bucket.

Multiplying three 66-limb difference magnitudes with the generic fixed-width
limb model naturally produces a 198-limb result type, providing 6336 bits.
That same width is sufficient for the final sign buckets.

The full-range fallback can therefore remain fixed-width and allocation-free.
This conservatively wide stack representation is acceptable because the
fallback is a correctness path rather than the ordinary hot path.

The intended binary64 execution path is therefore:

1. certified fast orient3d filter;
2. exact/adaptive expansion path for uncertain ordinary-range cases;
3. fixed-width exact dyadic fallback for extreme finite inputs.

The full-range dyadic path exists for correctness, not as the ordinary hot
path.

## Exact integral prototype evidence

The proposed fixed-width integral backend has been validated by a dedicated
research prototype under:

    experiments/orientation3_exact

The prototype is deliberately outside the production source tree.

It implements:

- exact signed coordinate differences as sign plus unsigned magnitude;
- 32-bit-limb fixed-width multiplication;
- positive and negative determinant-sign buckets;
- a 4-limb bucket for Point3!int;
- a 7-limb bucket for Point3!long;
- an independent std.bigint.BigInt determinant oracle.

### Explicit cases

The prototype verifies:

- the canonical positive determinant convention;
- exactly coplanar points;
- repeated-point degeneracy;
- complete int coordinate-span orientation;
- complete long coordinate-span orientation;
- exact coplanarity across complete coordinate spans;
- cancellation-heavy long coplanarity with large non-zero components.

### Permutation properties

Every checked geometry also verifies:

- one transposition reverses sign;
- a four-cycle reverses sign;
- two transpositions preserve sign.

### Exhaustive Leibniz-sign verification

All 2^9 = 512 assignments of non-zero signs to the nine components of the
three difference vectors were checked exhaustively.

The only possible positive/negative Leibniz-term splits are:

- 5 / 1;
- 3 / 3;
- 1 / 5.

The observed complete distribution is:

- 96 assignments with 5 / 1;
- 320 assignments with 3 / 3;
- 96 assignments with 1 / 5.

No sign bucket can therefore contain all six non-zero determinant terms.

This experimentally confirms the algebraic sign invariant used by the
fixed-width bucket bounds.

### Random exact-oracle comparison

For each compiler run, the optimized fixed-width predicate was compared
against the independent BigInt oracle for:

- 50,000 random Point3!int quadruples;
- 50,000 random Point3!long quadruples.

All comparisons passed.

### Compiler and optimization result

The complete integral prototype passes under:

- DMD debug;
- LDC debug;
- DMD release;
- LDC release.

Every configuration passes:

- the explicit semantic and extreme-value cases;
- the exhaustive 512-case Leibniz-sign test;
- the permutation properties;
- 50,000 random Point3!int oracle comparisons;
- 50,000 random Point3!long oracle comparisons.

Thus each compiler/configuration run performs 100,000 random exact-oracle
comparisons in addition to the deterministic cases.

### Research conclusion for integral coordinates

The exact integral architecture is considered validated for production design:

- Point3!int may use four 32-bit limbs for determinant sign buckets;
- Point3!long may use seven 32-bit limbs;
- signed arbitrary-width arithmetic is unnecessary;
- production BigInt is unnecessary;
- BigInt should remain an independent unittest oracle.

Any later production implementation should preserve these semantics but need
not copy the research prototype mechanically.

## Floating-point policy

### float

Every finite binary32 coordinate is exactly representable as binary64.

Point3!float may therefore promote exactly to Point3!double and use the
binary64 robust backend without loss of predicate information.

### double

Point3!double requires a certified robust predicate.

The intended architecture is staged:

1. fast ordinary binary64 determinant evaluation;
2. certified error bound;
3. return immediately when the sign is proven;
4. adaptive exact expansion fallback for uncertain cases;
5. exact full-range fallback for finite values not safely covered by the
   ordinary expansion path.

The architecture should follow the same correctness model as geo-d's robust
Orientation2 implementation while using a dimension-specific orient3d filter
and fallback.

A classical Shewchuk orient3d implementation is a primary algorithmic
reference, but its sign convention must be adapted to the geo3-d definition.

The final public guarantee is determinant-sign correctness, not closeness of a
computed determinant magnitude.

### Certified binary64 fast-filter evidence

The proposed first-stage binary64 orient3d filter has been validated by a
dedicated research prototype under:

    experiments/orientation3_double_filter

The filter uses:

- explicit binary64 rounding points;
- the classical orient3d difference layout relative to d;
- the Shewchuk first-stage orient3d error bound;
- conservative rejection of non-finite intermediates;
- conservative rejection of subnormal arithmetic;
- explicit detection of potentially underflowed products;
- sign normalization to geo3-d's public convention
  det(b-a, c-a, d-a).

The prototype compares every certified filter result against an independent
exact binary64-to-BigInt dyadic oracle.

The central validation rule is:

    uncertain is always permitted;
    every non-uncertain result must exactly match the oracle.

#### Ordinary-range random input

For 50,000 deterministic random moderate-coordinate quadruples:

    positive  = 24,888
    negative  = 25,112
    coplanar  = 0
    uncertain = 0

Thus every ordinary-range random case was certified by the fast path.

Every certified result matched the exact oracle.

#### Arbitrary finite binary64 input

For 5,000 deterministic quadruples generated from arbitrary finite binary64
bit patterns:

    positive  = 200
    negative  = 208
    coplanar  = 0
    uncertain = 4,592

The high uncertain rate is intentional.

Arbitrary binary64 bit patterns frequently produce exponent imbalance,
overflowing coordinate differences or products, subnormal intermediates, or
other cases outside the deliberately narrow first-stage proof.

The filter remains conservative rather than attempting to certify such cases.

Every one of the 408 certified results matched the exact oracle.

#### Exact coplanarity

For 10,000 random moderate-coordinate quadruples constrained to z = 0:

    positive  = 0
    negative  = 0
    coplanar  = 10,000
    uncertain = 0

Thus structurally exact zero determinants are recognized directly where no
unsafe arithmetic is involved.

#### Explicit extreme cases

The prototype also verifies:

- the canonical geo3-d positive sign convention;
- sign reversal under a transposition;
- a smallest-subnormal perturbation is returned as uncertain;
- overflowing determinant products are returned as uncertain;
- overflowing finite coordinate subtraction is returned as uncertain;
- NaN is never certified;
- infinity is never certified.

#### Compiler and optimization result

The complete filter experiment produces identical results under:

- DMD debug;
- LDC debug;
- DMD release;
- LDC release.

All four runs pass every oracle assertion.

### Research conclusion for the first binary64 stage

The certified first-stage filter design is considered validated.

Production design may use the same architecture:

1. reject unsupported/non-finite public input according to the public
   predicate contract;
2. evaluate a certified binary64 orient3d filter;
3. normalize the backend sign to det(b-a, c-a, d-a);
4. return immediately when the sign is certified;
5. route uncertain finite cases to an exact fallback.

The fast filter is not itself the complete Point3!double backend.

A robust public orientation operation must not expose the uncertain state.

### real

Point3!real orientation should remain unsupported initially.

It must not be implemented by demoting real to double.

Support may be added later only when a robust backend exists for the actual
platform representation of D real.

## Exact binary64 expansion evidence

The second binary64 stage has been validated by the research prototype:

    experiments/orientation3_double_expansion

This prototype evaluates the complete determinant exactly with floating-point
expansions inside a deliberately conservative exponent working range.

It uses exact snapshots of the numerical primitives from geo-d v2.0.0:

- repository: geo-d;
- tag: v2.0.0;
- commit: 83c974e0ca018bb378ce65a08ef9d3178decf9ad;
- expansion.d blob:
  50052c35af86e2d4306bf0e888cb5b37bb562c69;
- binary64_rounding.d blob:
  51c760f8fd4b7a58f6932ccfc16c635956aebcad.

The snapshots are research references only. They do not establish a production
dependency between geo3-d and geo-d.

### Expansion architecture

Each coordinate difference is represented exactly by a TwoDiff expansion of
at most two components.

The prototype then computes the six Leibniz triple products exactly:

    + ux * vy * wz
    + uy * vz * wx
    + uz * vx * wy
    - uz * vy * wx
    - uy * vx * wz
    - ux * vz * wy

Each triple product is constructed through expansion scaling and exact
expansion summation.

The six terms are then combined with a balanced exact expansion-sum tree and
the final orientation is obtained from the exact expansion sign.

No determinant magnitude is exposed as public API.

### Expansion test evidence

The prototype verifies explicitly:

- cancellation-heavy exact coplanarity;
- one-ULP perturbations on both sides of a plane;
- a non-zero TwoDiff tail;
- complete degeneracy;
- deliberate rejection of full-range overflowing differences;
- deliberate rejection of subnormal-scale cases outside the conservative
  expansion working range.

It additionally checks 50,000 deterministic moderate-coordinate random
quadruples against an independent exact BigInt dyadic oracle.

Every tested result matches the oracle.

### Expansion compiler and optimization result

The complete expansion prototype passes under:

- DMD debug;
- LDC debug;
- DMD release;
- LDC release.

All four runs produce:

    orientation3 exact expansion prototype:
    PASS (50000 moderate random cases)

The expansion stage is therefore considered validated for the conservative
ordinary exponent range.

Inputs outside its proven working range must fall through to the exact
full-range dyadic backend.

## Full-range binary64 dyadic evidence

The final binary64 stage has been validated by the research prototype:

    experiments/orientation3_double_dyadic

This stage supports the complete finite IEEE binary64 coordinate domain.

It uses exact snapshots from geo-d v2.0.0:

- fixed_uint.d blob:
  d50a7dc640584e3a4a95d6c49efa2007616b8104;
- dyadic.d blob:
  554a657f271da803dd6e524117d1e5bc6986090a.

Again, these are research references rather than a production dependency.

### Dyadic architecture

Every finite binary64 coordinate is decoded exactly as an integer multiple of:

    2^-1074

One coordinate magnitude requires 66 32-bit limbs.

Exact coordinate subtraction remains within the same 66-limb representation.

A product of two coordinate differences uses:

    132 limbs

and the 3D triple product uses:

    198 limbs
    = 6336 bits

The previously derived determinant bucket bound is less than 6300 bits.

Therefore the 198-limb representation is sufficient for:

- every exact triple product;
- every complete positive determinant bucket;
- every complete negative determinant bucket.

The determinant sign is obtained by comparing the two exact unsigned buckets.

No floating-point arithmetic occurs after coordinate decoding.

### Full-range explicit cases

The prototype verifies:

- the canonical positive public sign convention;
- exact coplanarity;
- complete degeneracy;
- the smallest positive binary64 subnormal in all three basis directions;
- determinants whose magnitude is far beyond binary64;
- finite coordinate subtraction that would overflow ordinary binary64;
- exact coplanarity across overflowing coordinate spans;
- extreme exponent imbalance;
- transposition sign reversal;
- four-cycle sign reversal.

### Full-range oracle evidence

Each run additionally checks:

- 5,000 arbitrary finite binary64 point quadruples;
- 5,000 arbitrary full-range quadruples constrained to z = 0.

Every general case is compared with an independent exact BigInt dyadic oracle.

Every constrained z = 0 case is also required explicitly to produce exact
coplanarity.

All tests pass.

### Full-range compiler and optimization result

The full-range dyadic prototype passes under:

- DMD debug;
- LDC debug;
- DMD release;
- LDC release.

All four configurations produce:

    orientation3 full-range dyadic prototype:
    PASS (5000 arbitrary finite + 5000 full-range coplanar cases)

The full finite binary64 domain is therefore covered by an experimentally
validated exact fallback.

## Validated numerical pipeline

The complete proposed Orientation3 numerical architecture is now covered by
independent research prototypes.

For integral coordinates:

    int
        -> exact fixed-width determinant sign

    long
        -> exact fixed-width determinant sign

For binary floating-point coordinates:

    float
        -> exact promotion to double
        -> binary64 robust pipeline

    double
        -> certified first-stage orient3d filter
        -> exact expansion fallback
        -> exact full-range dyadic fallback

The public predicate never exposes the internal uncertain state.

For every finite supported input, the final returned state is one of:

- negative;
- coplanar;
- positive.

Point3!real remains intentionally unsupported by the robust orientation
predicate.

## Non-finite floating input

As in geo-d Orientation2, robust floating orientation operates on finite
coordinates.

NaN and infinity are outside the mathematical predicate domain.

The public contract is therefore a finite-coordinate precondition rather
than a try-style result, matching the established Orientation2 family policy.

## No epsilon

Coplanarity must not be defined as

    abs(D) <= epsilon

for any global epsilon.

Near-coplanar but non-coplanar inputs must retain their mathematically correct
sign.

Applications that require a tolerance-based notion of approximate planarity
need a separate metric operation and must not alter Orientation3 semantics.

## Required independent tests

The implementation must be checked against independent exact oracles.

### Basic semantics

Test:

- canonical positive tetrahedron;
- canonical negative permutation;
- exactly coplanar points;
- repeated points;
- collinear base triples;
- fully degenerate four-point inputs.

### Permutations

For representative non-coplanar inputs:

- every transposition reverses orientation;
- every even permutation preserves orientation;
- every odd permutation reverses orientation.

### Integral extremes

Test complete-domain patterns around:

- int.min and int.max;
- long.min and long.max;
- exact coplanarity near full coordinate span;
- determinants requiring more than native 64 or 128-bit arithmetic.

Compare optimized predicates with an independent BigInt oracle.

### Floating near-degeneracy

Test:

- exactly coplanar binary64 inputs;
- the nearest representable perturbations away from coplanarity;
- subnormal coordinates;
- values near double.max;
- mixed very large and very small finite coordinates;
- cases forcing each adaptive fallback stage.

The final expected result should come from an independent exact dyadic or
BigInt-based test oracle rather than from ordinary floating arithmetic.

### Cross compiler

The predicate suite must pass at least:

- the supported minimum DMD;
- latest DMD;
- latest LDC.

Release builds and the existing lifetime gate remain required.

## Relationship to later 3D segment intersection

Orientation3 can establish coplanarity of four segment endpoints, but it does
not by itself solve three-dimensional segment intersection.

Later segment-intersection research must still handle:

- skew segments;
- collinear overlap;
- endpoint contact;
- degenerate segments;
- construction of an intersection point;
- representability and rounding of constructed coordinates.

Orientation3 should therefore be implemented as an independent foundational
predicate rather than designed around one intersection algorithm.

## Relationship to LinearRing3View

LinearRing3View remains a representation of a cyclic 3D point sequence and
does not require planarity.

Orientation3 may later help a planarity or topology-validation algorithm, but
the existence of Orientation3 does not make planarity an invariant of
LinearRing3View.

Ring-validation policy remains a separate research decision.

## Public API restraint

This research does not justify adding:

- Plane3;
- Triangle3;
- Tetrahedron3;
- Polygon3View;
- a generic public determinant API;
- a public cross-product API solely for orientation;
- generic public N-dimensional predicates;
- a tolerance parameter on orientation;
- approximate or fast public orientation variants.

Such API requires separate consumers or research.

## Research decision

The Orientation3 and coplanarity research is considered complete enough to
proceed to production implementation.

The public semantic contract should be:

    Orientation3 orientation(
        Point3!T a,
        Point3!T b,
        Point3!T c,
        Point3!T d
    )

for the supported robust scalar types.

The defining determinant is:

    det(b-a, c-a, d-a)

with:

    determinant < 0 -> Orientation3.negative
    determinant = 0 -> Orientation3.coplanar
    determinant > 0 -> Orientation3.positive

The result enum should be declared so that:

    Orientation3.init == Orientation3.coplanar

No public three-point Point3 orientation overload should be added.

No public coplanar(a, b, c, d) convenience function should be added initially.

No epsilon or tolerance participates in Orientation3 semantics.

Degenerate affine configurations are valid inputs and naturally produce
Orientation3.coplanar where the determinant is exactly zero.

## Supported scalar policy

The initial robust public overload set should support:

- Point3!int;
- Point3!long;
- Point3!float;
- Point3!double.

Point3!real should remain unsupported.

The scalar-specific numerical policy is:

### int

Use exact fixed-width integer arithmetic.

The validated determinant bucket width is four 32-bit limbs.

### long

Use exact fixed-width integer arithmetic.

The validated determinant bucket width is seven 32-bit limbs.

### float

Promote exactly to binary64 and use the double backend.

### double

Use the validated staged robust pipeline:

1. certified floating orient3d filter;
2. exact expansion fallback inside its proven working range;
3. exact fixed-width dyadic fallback for the remaining finite inputs.

Non-finite binary floating-point coordinates remain outside the public robust
predicate domain.

## Production implementation plan

Production implementation should proceed in small reviewable slices.

### Slice 1: public type and integral predicate

Add:

    source/geo3/orientation.d

containing:

- Orientation3;
- Point3!int orientation overload;
- Point3!long orientation overload.

Add private implementation support under:

    source/geo3/internal/

for the fixed-width exact determinant machinery.

Required tests include:

- canonical sign convention;
- coplanarity;
- all documented degeneracies;
- permutation identities;
- complete int and long coordinate extremes;
- independent BigInt oracle comparisons.

Do not expose the determinant magnitude.

### Slice 2: certified binary64 filter

Add the private orient3d first-stage filter.

Requirements:

- explicit binary64 rounding points;
- conservative uncertainty handling;
- exact sign normalization to the public determinant convention;
- no public uncertain state;
- filter decisions checked against an independent exact oracle.

The binary64 rounding helper may initially remain geo3-d-private.

Sharing implementation with geo-d must be considered separately from public
declaration sharing.

### Slice 3: exact expansion fallback

Add a geo3-d-private exact expansion fallback.

The implementation may be derived from the already validated expansion
architecture, but production source should not depend on geo-d internal
modules.

Before copying or adapting implementation code, preserve applicable repository
licensing and provenance information.

The production backend may later optimize the research implementation, but
must preserve its exact-sign contract.

### Slice 4: full-range dyadic fallback

Add the exact finite-binary64 fallback using:

- 66-limb exact coordinates/differences;
- 132-limb pair products;
- 198-limb triple products and determinant buckets.

The implementation must remain:

- fixed-width;
- allocation-free;
- pure;
- nothrow;
- @safe;
- @nogc.

A production BigInt dependency is not required.

### Slice 5: float overload

Add Point3!float by exact promotion to Point3!double.

No separate approximate float algorithm is required.

### Slice 6: package export and family audit

Export:

- Orientation3;
- orientation.

Then update:

    docs/research/api-family-audit.md

from research-required to aligned for robust 3D orientation/coplanarity.

The audit should continue to leave:

- 3D segment relationships/intersection;
- 3D ring validation;
- higher-level plane/surface abstractions

as separately gated research topics.

## Production quality gates

Before merging the production feature, require at least:

- git diff --check;
- DMD tests;
- LDC tests;
- DMD release build;
- LDC release build;
- existing DIP1000 lifetime gate;
- independent exact-oracle tests;
- extreme integer tests;
- near-coplanar binary64 tests;
- subnormal binary64 tests;
- complete finite-range binary64 fallback tests;
- permutation-property tests.

The research experiments should remain available until production tests cover
their important invariants.

## Architectural conclusion

Orientation3 is a genuine 3D counterpart to Orientation2 at the API-family
level, but its numerical implementation is substantially more demanding.

The established design rule remains:

    symmetry where the mathematics is symmetric;
    specialization where it is not.

The public family symmetry is therefore:

    orientation(a, b, c)       -> Orientation2
    orientation(a, b, c, d)    -> Orientation3

while the internal robust arithmetic remains dimension-specific.

No generic public N-dimensional orientation abstraction is justified by this
research.

## References reviewed

Primary references reviewed during this research:

- geo-d v2.0.0 source/geo/orientation.d
- geo-d ADR-0004, Numerical robustness and geometric predicates
- geo-d ADR-0016, Robust real scalar policy
- geo-d ADR-0020, Shared Euclidean contract core
- CGAL Kernel Orientation_3 documentation
- CGAL CoplanarOrientation_3 documentation
- Jonathan Richard Shewchuk, Adaptive Precision Floating-Point Arithmetic and
  Fast Robust Geometric Predicates
- Shewchuk public-domain predicates implementation and orient3d formulation
