module geo.internal.binary64_rounding;


/*
 * INTERNAL IMPLEMENTATION MODULE.
 *
 * Explicit binary64 rounding operations for robust floating-point
 * predicates.
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
 */


version (LDC)
{
    import ldc.llvmasm : __ir_pure;
}
else
{
    import core.math : toPrec;
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
        return toPrec!double(lhs + rhs);
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
        return toPrec!double(lhs - rhs);
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
        return toPrec!double(lhs * rhs);
    }
}


@safe unittest
{
    assert(roundedAdd(1.0, 2.0) == 3.0);
    assert(roundedSub(3.0, 2.0) == 1.0);
    assert(roundedMul(3.0, 2.0) == 6.0);

    assert(
        roundedAdd(
            1.0,
            0x1p-53
        ) == 1.0
    );

    assert(
        roundedSub(
            1.0,
            0x1p-53
        ) == 0x1.fffffffffffffp-1
    );
}
