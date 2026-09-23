/**
 * Explicit binary64 rounding operations for robust floating-point predicates.
 *
 * INTERNAL IMPLEMENTATION MODULE.
 *
 * Robust error bounds and error-free transformations depend on known
 * rounding points after individual elementary operations.
 *
 * LDC:
 *     Use plain LLVM binary64 fadd/fsub/fmul through inline IR.
 *     No fast-math flags are attached.
 *
 * Other compilers:
 *     Fall back to core.math.toPrec!double.
 *
 * Provenance:
 *     Derived from geo-d v2.0.0:
 *
 *         repository: https://github.com/alex-1974/geo-d
 *         tag:        v2.0.0
 *         commit:     83c974e0ca018bb378ce65a08ef9d3178decf9ad
 *         source:     source/geo/internal/binary64_rounding.d
 *         blob:       51c760f8fd4b7a58f6932ccfc16c635956aebcad
 *
 *     The implementation remains geo3-d-private; this does not create a
 *     production dependency on geo-d.
 *
 * Authors:
 *     Alexander Bernardi
 *
 * Copyright:
 *     Copyright © 2026 Alexander Bernardi
 *
 * License:
 *     MIT
 */
module geo3.internal.binary64_rounding;


version (LDC)
{
    import ldc.llvmasm :
        __ir_pure;
}
else
{
    import core.math :
        toPrec;
}


/**
 * Binary64 addition rounded at this operation.
 */
double roundedAdd(double lhs, double rhs)
    pure nothrow @safe @nogc
{
    version (LDC)
    {
        return __ir_pure!(
            `%r = fadd double %0, %1
             ret double %r`,
            double
        )(lhs, rhs);
    }
    else
    {
        return toPrec!double(
            lhs + rhs
        );
    }
}


/**
 * Binary64 subtraction rounded at this operation.
 */
double roundedSub(double lhs, double rhs)
    pure nothrow @safe @nogc
{
    version (LDC)
    {
        return __ir_pure!(
            `%r = fsub double %0, %1
             ret double %r`,
            double
        )(lhs, rhs);
    }
    else
    {
        return toPrec!double(
            lhs - rhs
        );
    }
}


/**
 * Binary64 multiplication rounded at this operation.
 */
double roundedMul(double lhs, double rhs)
    pure nothrow @safe @nogc
{
    version (LDC)
    {
        return __ir_pure!(
            `%r = fmul double %0, %1
             ret double %r`,
            double
        )(lhs, rhs);
    }
    else
    {
        return toPrec!double(
            lhs * rhs
        );
    }
}


@safe unittest
{
    assert(
        roundedAdd(
            1.0,
            2.0
        ) == 3.0
    );

    assert(
        roundedSub(
            3.0,
            2.0
        ) == 1.0
    );

    assert(
        roundedMul(
            3.0,
            2.0
        ) == 6.0
    );


    /*
     * Tie-to-even at an explicit binary64 rounding point.
     */
    assert(
        roundedAdd(
            1.0,
            0x1p-53
        ) == 1.0
    );


    /*
     * Adjacent representable value below 1.0.
     */
    assert(
        roundedSub(
            1.0,
            0x1p-53
        ) ==
        0x1.fffffffffffffp-1
    );
}
