# Orientation3 binary64 performance research

This experiment records the performance investigation that led to the
binary64 `Orientation3` first-stage filter optimization.

The production change is intentionally separate from this benchmark. The
benchmark exists to preserve the evidence and methodology used to evaluate
the optimization.

## Scope

The benchmark compares:

- the original production `Orientation3` pipeline;
- the original certified binary64 filter;
- an unchecked moderate-range reference filter;
- IEEE-754 bit-classification filter variants;
- an end-to-end candidate using the production exact fallbacks;
- the optimized production finite filter;
- an exact benchmark-local copy of the production filter, generated from the
  current production source before the benchmark is built.

The exact expansion and full-range dyadic fallback implementations are not
reimplemented here. Candidate pipelines use the production implementations.

## Correctness checks

The experiment includes:

- moderate-range equivalence checks;
- filter equivalence checks;
- full-range finite binary64 randomized checks;
- end-to-end candidate versus public `orientation` checks;
- near-degenerate cases;
- extreme-scale cases.

The optimized first-stage filter remains conservative: cases that cannot be
certified fall through to the exact expansion or dyadic backend.

## Main result

Replacing repeated floating-point classification helpers with direct
IEEE-754 binary64 classification substantially reduces the normal finite
`double` hot-path cost while preserving certified fallback semantics.

The production patch is:

    77c619e perf: optimize binary64 Orientation3 filter

No compiler-specific `Orientation3` inline policy is required by the
production implementation.

## LDC static-library boundary

A separate experiment compared:

1. `orientationFilterFinite` from the normally built `geo3-d` static library;
2. an exact benchmark-local copy of the same production source;
3. the experimental bit-guard filter.

The benchmark-local exact production copy matched the experimental filter,
while the separately compiled static-library function was slower.

`source/production_filter_local.d` is generated automatically from
`source/geo3/internal/orientation_filter.d` by `generate_local_filter.py`.
Only the module declaration is changed. The generated file is deliberately
not version-controlled so that the boundary experiment cannot silently drift
away from the current production implementation.

With:

    dub build --combined --build=release --compiler=ldc2

the production filter, its exact local copy, and the experimental filter
converged to essentially the same performance.

This isolates the remaining ordinary LDC build overhead to the compilation /
static-library boundary rather than to the `Orientation3` algorithm.

The general question of DUB combined builds and LDC LTO belongs to wider
D-language build-performance research and is intentionally not encoded as an
`Orientation3` compiler-specific optimization.

## Running

DMD:

    dub build --root=. --build=release --compiler=dmd --force
    taskset -c 2 ./orientation3-benchmark

LDC normal library build:

    dub build --root=. --build=release --compiler=ldc2 --force
    taskset -c 2 ./orientation3-benchmark

LDC combined build:

    dub build --root=. --build=release --compiler=ldc2 --combined --force
    taskset -c 2 ./orientation3-benchmark

Performance results are sensitive to CPU frequency, thermal state, and
benchmark ordering. Repeated runs and counterbalanced comparisons should be
used for conclusions.
