/**
 * Fixed-capacity binary64 expansion arithmetic for robust predicates.
 *
 * INTERNAL IMPLEMENTATION MODULE.
 *
 * This is deliberately not a general arbitrary-precision arithmetic
 * facility. It contains only the error-free transformations required by
 * geo3-d's robust Orientation3 expansion backend.
 *
 * Provenance:
 *     Derived from geo-d v2.0.0:
 *
 *         repository: https://github.com/alex-1974/geo-d
 *         tag:        v2.0.0
 *         commit:     83c974e0ca018bb378ce65a08ef9d3178decf9ad
 *         source:     source/geo/internal/expansion.d
 *         blob:       50052c35af86e2d4306bf0e888cb5b37bb562c69
 *
 *     The implementation remains geo3-d-private; this does not establish
 *     a production dependency on geo-d.
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
module geo3.internal.expansion;

import geo3.internal.binary64_rounding :
    roundedAdd,
    roundedMul,
    roundedSub;

import std.math.traits :
    isFinite;


/**
 * Exact two-component representation of one arithmetic result.
 *
 * Mathematically:
 *
 *     exact result = high + low
 */
struct TwoComponent
{
    double high;
    double low;
}


/*
 * High/low split of one binary64 value.
 */
private struct SplitComponent
{
    double high;
    double low;
}


/*
 * IEEE binary64 has 53 significant bits.
 *
 * Dekker / Shewchuk splitter:
 *
 *     2^ceil(53 / 2) + 1
 *     = 2^27 + 1
 */
private enum double splitter =
    134_217_729.0;


/**
 * Fixed-capacity inline storage for a floating-point expansion.
 *
 * Components are ordered from least significant to most significant.
 *
 * No allocation is performed.
 */
struct ExpansionBuffer(size_t Capacity)
if (Capacity > 0)
{
private:
    double[Capacity] _data;
    size_t _length;

public:
    enum size_t capacity =
        Capacity;


    @property size_t length() const
        pure nothrow @safe @nogc
    {
        return _length;
    }


    @property bool empty() const
        pure nothrow @safe @nogc
    {
        return _length == 0;
    }


    void clear()
        pure nothrow @safe @nogc
    {
        _length = 0;
    }


    void append(double value)
        pure nothrow @safe @nogc
    {
        assert(isFinite(value));
        assert(_length < Capacity);

        _data[_length] =
            value;

        ++_length;
    }


    double opIndex(size_t index) const
        pure nothrow @safe @nogc
    {
        assert(index < _length);

        return _data[index];
    }
}


/*
 * Absolute magnitude for a known finite value.
 */
private double finiteMagnitude(double value)
    pure nothrow @safe @nogc
{
    assert(isFinite(value));

    return value < 0.0
        ? -value
        : value;
}


/*
 * Error-free transformation of a + b.
 */
private TwoComponent twoSum(
    double a,
    double b
)
    pure nothrow @safe @nogc
{
    assert(isFinite(a));
    assert(isFinite(b));

    const double x =
        roundedAdd(
            a,
            b
        );

    assert(isFinite(x));

    const double bVirtual =
        roundedSub(
            x,
            a
        );

    const double aVirtual =
        roundedSub(
            x,
            bVirtual
        );

    const double bRoundoff =
        roundedSub(
            b,
            bVirtual
        );

    const double aRoundoff =
        roundedSub(
            a,
            aVirtual
        );

    const double y =
        roundedAdd(
            aRoundoff,
            bRoundoff
        );

    return TwoComponent(
        x,
        y
    );
}


/**
 * Error-free transformation of a - b.
 *
 * Preconditions:
 *
 * - inputs are finite;
 * - the rounded subtraction remains finite.
 *
 * Returns:
 *
 *     high + low == a - b
 *
 * exactly.
 */
TwoComponent twoDiff(
    double a,
    double b
)
    pure nothrow @safe @nogc
{
    assert(isFinite(a));
    assert(isFinite(b));

    const double x =
        roundedSub(
            a,
            b
        );

    assert(isFinite(x));

    const double bVirtual =
        roundedSub(
            a,
            x
        );

    const double aVirtual =
        roundedAdd(
            x,
            bVirtual
        );

    const double bRoundoff =
        roundedSub(
            bVirtual,
            b
        );

    const double aRoundoff =
        roundedSub(
            a,
            aVirtual
        );

    const double y =
        roundedAdd(
            aRoundoff,
            bRoundoff
        );

    return TwoComponent(
        x,
        y
    );
}


/*
 * Splits a binary64 value into non-overlapping high and low parts.
 *
 * The higher-level predicate is responsible for restricting inputs so that
 * multiplication by `splitter` cannot overflow or underflow critically.
 */
private SplitComponent split(double value)
    pure nothrow @safe @nogc
{
    assert(isFinite(value));

    const double c =
        roundedMul(
            splitter,
            value
        );

    assert(isFinite(c));

    const double aBig =
        roundedSub(
            c,
            value
        );

    const double high =
        roundedSub(
            c,
            aBig
        );

    const double low =
        roundedSub(
            value,
            high
        );

    return SplitComponent(
        high,
        low
    );
}


/*
 * Error-free transformation of a + b when |a| >= |b|.
 */
private TwoComponent fastTwoSum(
    double a,
    double b
)
    pure nothrow @safe @nogc
{
    assert(isFinite(a));
    assert(isFinite(b));

    assert(
        finiteMagnitude(a) >=
        finiteMagnitude(b)
    );

    const double x =
        roundedAdd(
            a,
            b
        );

    assert(isFinite(x));

    const double bVirtual =
        roundedSub(
            x,
            a
        );

    const double y =
        roundedSub(
            b,
            bVirtual
        );

    return TwoComponent(
        x,
        y
    );
}


/*
 * TwoProduct variant that reuses a precomputed split of b.
 */
private TwoComponent twoProductPresplit(
    double a,
    double b,
    SplitComponent bSplit
)
    pure nothrow @safe @nogc
{
    assert(isFinite(a));
    assert(isFinite(b));

    const double x =
        roundedMul(
            a,
            b
        );

    assert(isFinite(x));

    const SplitComponent aSplit =
        split(a);

    const double highProduct =
        roundedMul(
            aSplit.high,
            bSplit.high
        );

    const double err1 =
        roundedSub(
            x,
            highProduct
        );

    const double lowHighProduct =
        roundedMul(
            aSplit.low,
            bSplit.high
        );

    const double err2 =
        roundedSub(
            err1,
            lowHighProduct
        );

    const double highLowProduct =
        roundedMul(
            aSplit.high,
            bSplit.low
        );

    const double err3 =
        roundedSub(
            err2,
            highLowProduct
        );

    const double lowProduct =
        roundedMul(
            aSplit.low,
            bSplit.low
        );

    const double y =
        roundedSub(
            lowProduct,
            err3
        );

    return TwoComponent(
        x,
        y
    );
}


/*
 * Copies active non-zero components.
 *
 * Exact zero is canonicalized as [0.0].
 */
private void copyExpansionZeroElim(
    size_t SourceCapacity,
    size_t ResultCapacity
)(
    ref const ExpansionBuffer!SourceCapacity source,
    ref ExpansionBuffer!ResultCapacity result
)
    pure nothrow @safe @nogc
if (ResultCapacity >= SourceCapacity)
{
    result.clear();

    foreach (index; 0 .. source.length)
    {
        const double component =
            source[index];

        if (component != 0.0)
            result.append(component);
    }

    if (result.empty)
        result.append(0.0);
}


/**
 * Exact sum of two non-overlapping expansions with zero elimination.
 *
 * Inputs and result are ordered least-significant first.
 *
 * Maximum result length:
 *
 *     e.length + f.length
 */
void fastExpansionSumZeroElim(
    size_t ECapacity,
    size_t FCapacity,
    size_t HCapacity
)(
    ref const ExpansionBuffer!ECapacity e,
    ref const ExpansionBuffer!FCapacity f,
    ref ExpansionBuffer!HCapacity h
)
    pure nothrow @safe @nogc
if (HCapacity >= ECapacity + FCapacity)
{
    if (e.empty)
    {
        copyExpansionZeroElim(
            f,
            h
        );

        return;
    }

    if (f.empty)
    {
        copyExpansionZeroElim(
            e,
            h
        );

        return;
    }


    h.clear();

    size_t eIndex = 0;
    size_t fIndex = 0;

    double eNow =
        e[eIndex];

    double fNow =
        f[fIndex];

    double q;


    if (
        finiteMagnitude(eNow) <=
        finiteMagnitude(fNow)
    )
    {
        q =
            eNow;

        ++eIndex;
    }
    else
    {
        q =
            fNow;

        ++fIndex;
    }


    if (
        eIndex < e.length &&
        fIndex < f.length
    )
    {
        double next;

        if (
            finiteMagnitude(e[eIndex]) <=
            finiteMagnitude(f[fIndex])
        )
        {
            next =
                e[eIndex];

            ++eIndex;
        }
        else
        {
            next =
                f[fIndex];

            ++fIndex;
        }


        const TwoComponent initial =
            fastTwoSum(
                next,
                q
            );

        if (initial.low != 0.0)
            h.append(initial.low);

        q =
            initial.high;


        while (
            eIndex < e.length &&
            fIndex < f.length
        )
        {
            if (
                finiteMagnitude(e[eIndex]) <=
                finiteMagnitude(f[fIndex])
            )
            {
                next =
                    e[eIndex];

                ++eIndex;
            }
            else
            {
                next =
                    f[fIndex];

                ++fIndex;
            }


            const TwoComponent sum =
                twoSum(
                    q,
                    next
                );

            if (sum.low != 0.0)
                h.append(sum.low);

            q =
                sum.high;
        }
    }


    while (eIndex < e.length)
    {
        const TwoComponent sum =
            twoSum(
                q,
                e[eIndex]
            );

        if (sum.low != 0.0)
            h.append(sum.low);

        q =
            sum.high;

        ++eIndex;
    }


    while (fIndex < f.length)
    {
        const TwoComponent sum =
            twoSum(
                q,
                f[fIndex]
            );

        if (sum.low != 0.0)
            h.append(sum.low);

        q =
            sum.high;

        ++fIndex;
    }


    if (
        q != 0.0 ||
        h.empty
    )
        h.append(q);
}


/**
 * Multiplies an expansion by one binary64 scalar exactly and eliminates
 * zero components.
 *
 * Maximum result length:
 *
 *     2 * expansion.length
 *
 * The higher-level predicate constrains the working exponent range so the
 * EFT preconditions remain valid.
 */
void scaleExpansionZeroElim(
    size_t ECapacity,
    size_t HCapacity
)(
    ref const ExpansionBuffer!ECapacity expansion,
    double scalar,
    ref ExpansionBuffer!HCapacity result
)
    pure nothrow @safe @nogc
if (HCapacity >= 2 * ECapacity)
{
    assert(isFinite(scalar));

    if (expansion.empty)
    {
        result.clear();
        result.append(0.0);

        return;
    }


    result.clear();

    const SplitComponent scalarSplit =
        split(scalar);

    TwoComponent product =
        twoProductPresplit(
            expansion[0],
            scalar,
            scalarSplit
        );

    double q =
        product.high;

    if (product.low != 0.0)
        result.append(product.low);


    foreach (index; 1 .. expansion.length)
    {
        product =
            twoProductPresplit(
                expansion[index],
                scalar,
                scalarSplit
            );

        const TwoComponent sum =
            twoSum(
                q,
                product.low
            );

        if (sum.low != 0.0)
            result.append(sum.low);


        const TwoComponent accumulated =
            fastTwoSum(
                product.high,
                sum.high
            );

        if (accumulated.low != 0.0)
            result.append(
                accumulated.low
            );

        q =
            accumulated.high;
    }


    if (
        q != 0.0 ||
        result.empty
    )
        result.append(q);
}


/**
 * Exact sign of a valid non-overlapping expansion.
 *
 * Returns:
 *
 *     -1 negative
 *      0 exact zero
 *      1 positive
 */
int expansionSign(size_t Capacity)(
    ref const ExpansionBuffer!Capacity expansion
)
    pure nothrow @safe @nogc
{
    size_t index =
        expansion.length;

    while (index != 0)
    {
        --index;

        const double component =
            expansion[index];

        if (component > 0.0)
            return 1;

        if (component < 0.0)
            return -1;
    }

    return 0;
}


/**
 * Exact component-wise expansion negation.
 */
void negateExpansion(
    size_t SourceCapacity,
    size_t ResultCapacity
)(
    ref const ExpansionBuffer!SourceCapacity source,
    ref ExpansionBuffer!ResultCapacity result
)
    pure nothrow @safe @nogc
if (ResultCapacity >= SourceCapacity)
{
    result.clear();

    foreach (index; 0 .. source.length)
    {
        result.append(
            -source[index]
        );
    }
}


@safe unittest
{
    /*
     * TwoDiff must preserve a low-order value lost by the rounded result.
     */
    const TwoComponent difference =
        twoDiff(
            10_000_000_000_000_000.0,
            1.0
        );

    assert(
        difference.high ==
        10_000_000_000_000_000.0
    );

    assert(
        difference.low ==
        -1.0
    );


    /*
     * Scaling exercises the exact product tail.
     *
     * (1 + 2^-52)^2
     *   = 1 + 2^-51 + 2^-104
     */
    enum double onePlusUlp =
        0x1.0000000000001p+0;

    ExpansionBuffer!1 source;
    ExpansionBuffer!2 scaled;

    source.append(
        onePlusUlp
    );

    scaleExpansionZeroElim(
        source,
        onePlusUlp,
        scaled
    );

    assert(
        scaled.length ==
        2
    );

    assert(
        scaled[0] ==
        0x1p-104
    );

    assert(
        scaled[1] ==
        0x1.0000000000002p+0
    );


    /*
     * Exact cancellation removes low components.
     */
    ExpansionBuffer!2 lhs;
    ExpansionBuffer!2 rhs;
    ExpansionBuffer!4 sum;

    lhs.append(
        0x1p-104
    );

    lhs.append(
        1.0
    );

    rhs.append(
        -0x1p-104
    );

    rhs.append(
        2.0
    );

    fastExpansionSumZeroElim(
        lhs,
        rhs,
        sum
    );

    assert(
        sum.length ==
        1
    );

    assert(
        sum[0] ==
        3.0
    );


    /*
     * Exact sign uses the most significant non-zero component.
     */
    ExpansionBuffer!4 positive;

    positive.append(
        -0x1p-104
    );

    positive.append(
        0x1p-52
    );

    positive.append(
        1.0
    );

    assert(
        expansionSign(positive) ==
        1
    );


    ExpansionBuffer!4 negative;

    negateExpansion(
        positive,
        negative
    );

    assert(
        expansionSign(negative) ==
        -1
    );


    static assert(
        TwoComponent.sizeof ==
        2 * double.sizeof
    );

    static assert(
        ExpansionBuffer!4.capacity ==
        4
    );
}
