# containers-d StaticVector consumer proof

Status: qualified research evidence  
Branch: `research/containers-static-vector-proof`  
Draft PR: #29  
containers-d research: issue #25 / PR #33

## Purpose

This is the independent second-consumer proof for the containers-d fixed-vector
family research.

geo3-d contains its own robust binary64 expansion implementation for the
Orientation3 exact fallback. Its original `ExpansionBuffer!Capacity` is an
independent duplicate of the same fixed inline storage mechanics used by geo-d.

The experiment asks whether the same containers-d family primitive can replace
that duplicate without introducing a 2D-specific assumption or a performance
regression.

## Research integration

geo3-d keeps its numerical domain type:

```d
struct ExpansionBuffer(size_t Capacity)
```

and composes the research-only containers-d scalar fixed-vector mechanics:

```d
mixin ScalarStaticVectorOps!(double, Capacity);
```

geo3-d retains the only domain-specific operation:

```d
void append(double value)
    pure nothrow @safe @nogc
{
    assert(isFinite(value));
    pushBack(value);
}
```

The containers-d research dependency is pinned to:

`f853f6cb2c77b2a5afeebfb9bb2be30768422f80`

No containers-d type leaks into geo3-d's public API.

## Functional qualification

The branch passes geo3-d's complete ordinary CI on both baseline compiler
families:

- DMD 2.111;
- LDC 1.41;
- unit tests;
- lifetime compile-negative tests;
- release builds.

The expansion and Orientation3 semantic checks in the branch-vs-develop probe
also produce identical checksums.

## Real branch-vs-develop performance probe

Source:

`experiments/static_vector_consumer_probe.d`

Workflow:

`.github/workflows/perf-static-vector-consumer.yml`

The workflow builds exactly the same probe source twice on the same runner:

1. against a fresh checkout of `geo3-d develop`;
2. against the research branch plus the exact containers-d candidate.

The probe measures 131,072 calls for:

- `scaleExpansionZeroElim` with a 2-component input and 4-component result;
- `fastExpansionSumZeroElim` with two Capacity-4 inputs and Capacity-8 result;
- exact coplanar Orientation3 through `tryOrientationExactExpansion`;
- near-coplanar Orientation3 through the same exact expansion backend.

The Orientation3 backend exercises much larger intermediate buffers, including
`ExpansionBuffer!32`, so this is a stronger capacity/algorithm test than the
initial Capacity-1..4 mechanics probes.

## DMD 2.111 result

| Operation | develop Ir | composition Ir | delta |
|---|---:|---:|---:|
| scaleExpansion 2→4 | 30,933,012 | 30,933,012 | 0 |
| fastExpansionSum 4+4 | 40,501,268 | 40,501,268 | 0 |
| Orientation3 exact coplanar | 598,081,570 | 598,081,570 | 0 |
| Orientation3 exact near | 597,360,674 | 597,360,674 | 0 |

Every measured DMD path is instruction-identical to develop.

This independently reproduces the geo-d result and includes the larger
Orientation3 expansion graph.

## LDC 1.41 result

| Operation | develop Ir | composition Ir | delta |
|---|---:|---:|---:|
| scaleExpansion 2→4 | 15,728,671 | 15,728,671 | 0 |
| fastExpansionSum 4+4 | 37,290,019 | 37,552,163 | +262,144 |
| Orientation3 exact coplanar | 477,495,322 | 477,495,322 | 0 |
| Orientation3 exact near | 476,708,890 | 476,708,890 | 0 |

The only non-zero delta is:

```text
262,144 / 131,072 calls = +2 Ir/op
```

for `fastExpansionSum 4+4`.

This is the same +2 Ir/op LDC code-generation effect observed independently in
geo-d's corresponding sum probe.

It is therefore reproducible compiler behaviour, but still too small to
justify a compiler-specific consumer policy without evidence from a material
end-to-end workload.

## Family-model conclusion

The two independent geometry consumers now demonstrate the same successful
shape:

```text
containers-d
    generic scalar fixed-vector mechanics
                |
                | typed compile-time composition
                v
    +---------------------------+
    |                           |
geo-d ExpansionBuffer       geo3-d ExpansionBuffer
    |                           |
finite-component invariant  finite-component invariant
2D robust predicates        3D robust predicates
```

Both domain libraries keep their own semantic type and numerical invariants.
Neither pays runtime abstraction or policy-dispatch cost.

On DMD 2.111 the real measured consumer paths are instruction-identical to
their hand-local develop implementations.

On LDC 1.41 all measured paths are identical except the same reproducible
+2 Ir/op expansion-sum code-shape effect.

## Decision

The independent second-consumer gate for M4.3 is passed.

This supports promoting the **architecture** of consumer-side typed composition
as a valid containers-d adaptation mechanism.

It does not yet make `ScalarStaticVectorOps` a stable public API. Naming,
scope, relationship to the ordinary `StaticVector!(T,N)` type, diagnostics,
and long-term support still belong to the M4.3/M4.6 API decision.
