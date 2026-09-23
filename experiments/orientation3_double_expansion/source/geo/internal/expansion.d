module geo.internal.expansion;

import geo.internal.binary64_rounding :
    roundedAdd,
    roundedMul,
    roundedSub;
import std.math.traits : isFinite;


/*
 * INTERNAL IMPLEMENTATION MODULE.
 *
 * Error-free floating-point transformations used by robust predicates.
 *
 * This first implementation deliberately targets binary64 only.
 * `real` requires a separate precision-aware backend.
 */


/**
 * Two-component exact representation of one arithmetic result.
 *
 * Mathematically:
 *
 *     exact result = high + low
 *
 * `high` is the ordinary rounded binary64 result and `low` is the
 * exactly recovered roundoff term, subject to the documented
 * preconditions of the corresponding transformation.
 */
struct TwoComponent
{
    double high;
    double low;
}


/**
 * High/low split of one binary64 value.
 *
 * Mathematically:
 *
 *     value = high + low
 */
struct SplitComponent
{
    double high;
    double low;
}


/*
 * For IEEE binary64:
 *
 *     p = 53 significant bits
 *
 * Shewchuk / Dekker splitter:
 *
 *     2^ceil(p / 2) + 1
 *     = 2^27 + 1
 *     = 134217729
 */
private enum double splitter =
    134_217_729.0;

/**
 * Error-free transformation of a + b.
 *
 * Returns high and low such that, mathematically:
 *
 *     high + low == a + b
 *
 * while `high` equals the ordinary correctly rounded binary64 sum.
 *
 * The current robust-predicate backend calls this only where the
 * intermediate binary64 operations remain finite.
 */
TwoComponent twoSum(double a, double b)
    pure nothrow @safe @nogc
{
    assert(isFinite(a));
    assert(isFinite(b));

    const double x =
        roundedAdd(a, b);

    assert(isFinite(x));

    const double bVirtual =
        roundedSub(x, a);

    const double aVirtual =
        roundedSub(x, bVirtual);

    const double bRoundoff =
        roundedSub(b, bVirtual);

    const double aRoundoff =
        roundedSub(a, aVirtual);

    const double y =
        roundedAdd(
            aRoundoff,
            bRoundoff
        );

    return TwoComponent(x, y);
}


/**
 * Error-free transformation of a - b.
 *
 * Returns high and low such that, mathematically:
 *
 *     high + low == a - b
 *
 * while `high` equals the ordinary correctly rounded binary64
 * difference.
 *
 * The current robust-predicate backend calls this only where the
 * intermediate binary64 operations remain finite.
 */
TwoComponent twoDiff(double a, double b)
    pure nothrow @safe @nogc
{
    assert(isFinite(a));
    assert(isFinite(b));

    const double x =
        roundedSub(a, b);

    assert(isFinite(x));

    const double bVirtual =
        roundedSub(a, x);

    const double aVirtual =
        roundedAdd(x, bVirtual);

    const double bRoundoff =
        roundedSub(bVirtual, b);

    const double aRoundoff =
        roundedSub(a, aVirtual);

    const double y =
        roundedAdd(
            aRoundoff,
            bRoundoff
        );

    return TwoComponent(x, y);
}


/**
 * Splits a binary64 value into non-overlapping high and low parts.
 *
 * The caller must ensure that multiplication by `splitter` does not
 * overflow. The later robust-predicate layer may scale coordinates
 * before entering expansion arithmetic when necessary.
 */
SplitComponent split(double value)
    pure nothrow @safe @nogc
{
    assert(isFinite(value));

    const double c =
        roundedMul(
            splitter,
            value
        );

    /*
     * This assertion makes the current numerical domain explicit.
     * It is not a general overflow-recovery mechanism.
     */
    assert(isFinite(c));

    const double aBig =
        roundedSub(c, value);

    const double high =
        roundedSub(c, aBig);

    const double low =
        roundedSub(value, high);

    return SplitComponent(high, low);
}


/**
 * Error-free transformation of a * b.
 *
 * Returns high and low such that, mathematically:
 *
 *     high + low == a * b
 *
 * while `high` equals the ordinary rounded binary64 product.
 *
 * Preconditions of this initial backend:
 *
 * - a and b are finite;
 * - a * b does not overflow;
 * - splitting either operand does not overflow.
 *
 * Exponent scaling for inputs outside this safe working range belongs
 * to the higher-level robust-predicate implementation.
 */
TwoComponent twoProduct(double a, double b)
    pure nothrow @safe @nogc
{
    assert(isFinite(a));
    assert(isFinite(b));

    const double x =
        roundedMul(a, b);

    assert(isFinite(x));

    const auto aSplit =
        split(a);

    const auto bSplit =
        split(b);

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

    return TwoComponent(x, y);
}


/**
 * Fixed-capacity stack/value storage for a floating-point expansion.
 *
 * Components are stored from least significant to most significant.
 *
 * This type deliberately owns its storage inline:
 *
 * - no dynamic array allocation;
 * - no GC dependency;
 * - capacity is known at compile time.
 *
 * The type does not itself enforce expansion non-overlap or ordering.
 * Those invariants are established by the algorithms that produce an
 * expansion.
 */
struct ExpansionBuffer(size_t Capacity)
if (Capacity > 0)
{
private:
    double[Capacity] _data;
    size_t _length;

public:
    enum size_t capacity = Capacity;


    /// Number of active expansion components.
    @property size_t length() const
        pure nothrow @safe @nogc
    {
        return _length;
    }


    /// True when the expansion contains no active components.
    @property bool empty() const
        pure nothrow @safe @nogc
    {
        return _length == 0;
    }


    /**
     * Removes all active components.
     *
     * Stored bytes need not be cleared because values beyond `length`
     * are not part of the expansion.
     */
    void clear()
        pure nothrow @safe @nogc
    {
        _length = 0;
    }


    /**
     * Appends one component.
     *
     * Preconditions:
     *
     * - value is finite;
     * - spare capacity is available.
     */
    void append(double value)
        pure nothrow @safe @nogc
    {
        assert(isFinite(value));
        assert(_length < Capacity);

        _data[_length] = value;
        ++_length;
    }


    /**
     * Indexed access to active components.
     */
    double opIndex(size_t index) const
        pure nothrow @safe @nogc
    {
        assert(index < _length);
        return _data[index];
    }
}


/*
 * Absolute value for already finite binary64 inputs.
 */
private double finiteMagnitude(double value)
    pure nothrow @safe @nogc
{
    assert(isFinite(value));

    return value < 0.0
        ? -value
        : value;
}


/**
 * Error-free transformation of a + b when |a| >= |b|.
 *
 * Compared with twoSum(), FastTwoSum needs fewer operations because its
 * magnitude precondition guarantees the required rounding relation.
 *
 * Returns:
 *
 *     high + low == a + b
 *
 * mathematically, while `high` is the ordinary rounded binary64 sum.
 */
TwoComponent fastTwoSum(double a, double b)
    pure nothrow @safe @nogc
{
    assert(isFinite(a));
    assert(isFinite(b));

    assert(
        finiteMagnitude(a) >=
        finiteMagnitude(b)
    );

    const double x =
        roundedAdd(a, b);

    assert(isFinite(x));

    const double bVirtual =
        roundedSub(x, a);

    const double y =
        roundedSub(b, bVirtual);

    return TwoComponent(x, y);
}


/**
 * Floating-point estimate of an expansion's numerical value.
 *
 * Expansion components are accumulated from least significant to most
 * significant.
 *
 * This function is intentionally approximate. It must not be used as a
 * replacement for an exact expansion sign test.
 */
double estimate(size_t Capacity)(
    ref const ExpansionBuffer!Capacity expansion
)
    pure nothrow @safe @nogc
{
    double result = 0.0;

    foreach (index; 0 .. expansion.length)
    {
        result = roundedAdd(
            result,
            expansion[index]
        );
    }

    return result;
}


/*
 * TwoProduct variant that reuses an already computed split of b.
 *
 * This is the primitive required by scaleExpansionZeroElim().
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
        roundedMul(a, b);

    assert(isFinite(x));

    const auto aSplit =
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

    return TwoComponent(x, y);
}


/*
 * Copies active non-zero components.
 *
 * A numerically zero expansion is canonicalized as one zero component.
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
 * Adds two non-overlapping expansions and eliminates zero components.
 *
 * Input components must be ordered from least significant to most
 * significant.
 *
 * The result is likewise ordered from least significant to most
 * significant and contains no zero components unless the exact result
 * is zero, in which case the canonical result is:
 *
 *     [0.0]
 *
 * The output buffer must not alias either input buffer.
 *
 * Maximum output length:
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
        copyExpansionZeroElim(f, h);
        return;
    }

    if (f.empty)
    {
        copyExpansionZeroElim(e, h);
        return;
    }

    h.clear();

    size_t eIndex = 0;
    size_t fIndex = 0;

    double eNow = e[eIndex];
    double fNow = f[fIndex];

    double q;

    /*
     * Merge the components by increasing magnitude.
     */
    if (finiteMagnitude(eNow) <=
        finiteMagnitude(fNow))
    {
        q = eNow;
        ++eIndex;
    }
    else
    {
        q = fNow;
        ++fIndex;
    }

    /*
     * The second component is guaranteed to have magnitude >= |q|,
     * hence FastTwoSum is valid.
     */
    if (eIndex < e.length &&
        fIndex < f.length)
    {
        double next;

        if (finiteMagnitude(e[eIndex]) <=
            finiteMagnitude(f[fIndex]))
        {
            next = e[eIndex];
            ++eIndex;
        }
        else
        {
            next = f[fIndex];
            ++fIndex;
        }

        const auto initial =
            fastTwoSum(next, q);

        if (initial.low != 0.0)
            h.append(initial.low);

        q = initial.high;

        /*
         * Once Q has accumulated several components, the general
         * TwoSum transform is required.
         */
        while (eIndex < e.length &&
               fIndex < f.length)
        {
            if (finiteMagnitude(e[eIndex]) <=
                finiteMagnitude(f[fIndex]))
            {
                next = e[eIndex];
                ++eIndex;
            }
            else
            {
                next = f[fIndex];
                ++fIndex;
            }

            const auto sum =
                twoSum(q, next);

            if (sum.low != 0.0)
                h.append(sum.low);

            q = sum.high;
        }
    }

    while (eIndex < e.length)
    {
        const auto sum =
            twoSum(
                q,
                e[eIndex]
            );

        if (sum.low != 0.0)
            h.append(sum.low);

        q = sum.high;
        ++eIndex;
    }

    while (fIndex < f.length)
    {
        const auto sum =
            twoSum(
                q,
                f[fIndex]
            );

        if (sum.low != 0.0)
            h.append(sum.low);

        q = sum.high;
        ++fIndex;
    }

    if (q != 0.0 || h.empty)
        h.append(q);
}


/**
 * Multiplies an expansion by one binary64 scalar and eliminates zero
 * components.
 *
 * Input components must be ordered from least significant to most
 * significant and form a valid non-overlapping expansion.
 *
 * The result preserves expansion order and contains no zero components
 * unless the exact result is zero, in which case the canonical result
 * is:
 *
 *     [0.0]
 *
 * The output buffer must not alias the input buffer.
 *
 * Maximum output length:
 *
 *     2 * expansion.length
 *
 * This initial implementation inherits the current split()/TwoProduct
 * working-range preconditions. Exponent scaling for extreme values
 * belongs to the later predicate layer.
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

    const auto scalarSplit =
        split(scalar);

    auto product =
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

        const auto sum =
            twoSum(
                q,
                product.low
            );

        if (sum.low != 0.0)
            result.append(sum.low);

        /*
         * Expansion ordering guarantees the magnitude relation required
         * by FastTwoSum here.
         */
        const auto accumulated =
            fastTwoSum(
                product.high,
                sum.high
            );

        if (accumulated.low != 0.0)
            result.append(accumulated.low);

        q = accumulated.high;
    }

    if (q != 0.0 || result.empty)
        result.append(q);
}


/**
 * Exact sign of a valid non-overlapping expansion.
 *
 * Components are stored from least significant to most significant.
 * Therefore the highest-index non-zero component determines the exact
 * sign of the represented value.
 *
 * Returns:
 *
 *     -1  negative
 *      0  exact zero
 *      1  positive
 *
 * Unlike estimate(), this function is suitable for robust predicate
 * decisions.
 */
int expansionSign(size_t Capacity)(
    ref const ExpansionBuffer!Capacity expansion
)
    pure nothrow @safe @nogc
{
    size_t index = expansion.length;

    while (index > 0)
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
 * Component-wise negation of an expansion.
 *
 * Expansion order and non-overlap are preserved because multiplication
 * by -1 is exact in binary floating point.
 *
 * The output buffer must not alias the input buffer.
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
     * Exact expansion sign is determined by the most significant
     * non-zero component, not by a rounded estimate.
     */
    {
        ExpansionBuffer!4 positive;

        positive.append(-0x1p-104);
        positive.append(0x1p-52);
        positive.append(1.0);

        assert(expansionSign(positive) == 1);

        ExpansionBuffer!4 negative;

        negative.append(0x1p-104);
        negative.append(-0x1p-52);
        negative.append(-1.0);

        assert(expansionSign(negative) == -1);
    }


    /*
     * Zero components are skipped defensively.
     */
    {
        ExpansionBuffer!4 expansion;

        expansion.append(1.0);
        expansion.append(0.0);
        expansion.append(0.0);

        assert(expansionSign(expansion) == 1);
    }


    /*
     * Empty and canonical-zero expansions both represent exact zero.
     */
    {
        ExpansionBuffer!2 emptyExpansion;

        assert(
            expansionSign(emptyExpansion) == 0
        );

        ExpansionBuffer!2 zero;

        zero.append(0.0);

        assert(
            expansionSign(zero) == 0
        );
    }


    /*
     * Negation preserves component order and reverses exact sign.
     */
    {
        ExpansionBuffer!3 source;
        ExpansionBuffer!3 negated;

        source.append(0x1p-104);
        source.append(-0x1p-52);
        source.append(2.0);

        negateExpansion(
            source,
            negated
        );

        assert(negated.length == 3);

        assert(
            negated[0] == -0x1p-104
        );

        assert(
            negated[1] == 0x1p-52
        );

        assert(
            negated[2] == -2.0
        );

        assert(
            expansionSign(source) == 1
        );

        assert(
            expansionSign(negated) == -1
        );
    }


    /*
     * Double negation reconstructs the original expansion exactly.
     */
    {
        ExpansionBuffer!3 source;
        ExpansionBuffer!3 first;
        ExpansionBuffer!3 second;

        source.append(-0x1p-104);
        source.append(0x1p-52);
        source.append(1.0);

        negateExpansion(
            source,
            first
        );

        negateExpansion(
            first,
            second
        );

        assert(
            second.length ==
            source.length
        );

        foreach (index; 0 .. source.length)
        {
            assert(
                second[index] ==
                source[index]
            );
        }
    }


    /*
     * Fast expansion sum: exact cancellation removes zero components.
     */
    {
        ExpansionBuffer!2 e;
        ExpansionBuffer!2 f;
        ExpansionBuffer!4 h;

        e.append(0x1p-104);
        e.append(1.0);

        f.append(-0x1p-104);
        f.append(2.0);

        fastExpansionSumZeroElim(
            e,
            f,
            h
        );

        assert(h.length == 1);
        assert(h[0] == 3.0);
        assert(estimate(h) == 3.0);
    }


    /*
     * Fast expansion sum preserves a genuine low-order component.
     *
     * Exact value:
     *
     *     1 + 2 + 2^-104 + 2^-103
     *   = 3 + 3 * 2^-104
     */
    {
        ExpansionBuffer!2 e;
        ExpansionBuffer!2 f;
        ExpansionBuffer!4 h;

        e.append(0x1p-104);
        e.append(1.0);

        f.append(0x1p-103);
        f.append(2.0);

        fastExpansionSumZeroElim(
            e,
            f,
            h
        );

        assert(h.length == 2);

        assert(
            h[0] ==
            0x1.8p-103
        );

        assert(h[1] == 3.0);
    }


    /*
     * Empty inputs are supported internally and canonicalize exact zero
     * to one zero component.
     */
    {
        ExpansionBuffer!1 e;
        ExpansionBuffer!1 f;
        ExpansionBuffer!2 h;

        fastExpansionSumZeroElim(
            e,
            f,
            h
        );

        assert(h.length == 1);
        assert(h[0] == 0.0);
    }


    /*
     * Adding an empty expansion copies the non-zero expansion.
     */
    {
        ExpansionBuffer!2 e;
        ExpansionBuffer!2 emptyExpansion;
        ExpansionBuffer!4 h;

        e.append(0x1p-100);
        e.append(1.0);

        fastExpansionSumZeroElim(
            e,
            emptyExpansion,
            h
        );

        assert(h.length == 2);
        assert(h[0] == 0x1p-100);
        assert(h[1] == 1.0);
    }


    /*
     * Scaling by an exact power of two preserves the expansion shape.
     */
    {
        ExpansionBuffer!2 e;
        ExpansionBuffer!4 h;

        e.append(0x1p-104);
        e.append(1.0);

        scaleExpansionZeroElim(
            e,
            2.0,
            h
        );

        assert(h.length == 2);
        assert(h[0] == 0x1p-103);
        assert(h[1] == 2.0);
    }


    /*
     * Scaling one component exercises the exact TwoProduct tail.
     *
     * Let:
     *
     *     a = 1 + 2^-52
     *
     * Then:
     *
     *     a² = 1 + 2^-51 + 2^-104
     */
    {
        enum double a =
            0x1.0000000000001p+0;

        ExpansionBuffer!1 e;
        ExpansionBuffer!2 h;

        e.append(a);

        scaleExpansionZeroElim(
            e,
            a,
            h
        );

        assert(h.length == 2);

        assert(
            h[0] ==
            0x1p-104
        );

        assert(
            h[1] ==
            0x1.0000000000002p+0
        );
    }


    /*
     * Scaling by zero produces the canonical zero expansion.
     */
    {
        ExpansionBuffer!2 e;
        ExpansionBuffer!4 h;

        e.append(0x1p-104);
        e.append(1.0);

        scaleExpansionZeroElim(
            e,
            0.0,
            h
        );

        assert(h.length == 1);
        assert(h[0] == 0.0);
    }


    /*
     * ExpansionBuffer owns fixed inline storage and begins empty.
     */
    {
        ExpansionBuffer!4 expansion;

        assert(expansion.empty);
        assert(expansion.length == 0);
        static assert(
            ExpansionBuffer!4.capacity == 4
        );

        assert(estimate(expansion) == 0.0);

        expansion.append(0x1p-104);
        expansion.append(0x1p-52);
        expansion.append(1.0);

        assert(!expansion.empty);
        assert(expansion.length == 3);

        assert(expansion[0] == 0x1p-104);
        assert(expansion[1] == 0x1p-52);
        assert(expansion[2] == 1.0);

        expansion.clear();

        assert(expansion.empty);
        assert(expansion.length == 0);
        assert(estimate(expansion) == 0.0);
    }


    /*
     * FastTwoSum recovers an addend lost from the rounded high
     * component.
     */
    {
        enum double halfUlpAtOne =
            0x1p-53;

        const auto r =
            fastTwoSum(
                1.0,
                halfUlpAtOne
            );

        assert(r.high == 1.0);
        assert(r.low == halfUlpAtOne);
    }


    /*
     * FastTwoSum also handles exact additions naturally.
     */
    {
        const auto r =
            fastTwoSum(4.0, 2.0);

        assert(r.high == 6.0);
        assert(r.low == 0.0);
    }


    /*
     * Sign handling for FastTwoSum.
     */
    {
        const auto r =
            fastTwoSum(
                -1.0,
                -0x1p-53
            );

        assert(r.high == -1.0);
        assert(r.low == -0x1p-53);
    }


    /*
     * Cancellation remains exact when the magnitude precondition holds.
     */
    {
        const auto r =
            fastTwoSum(
                1.0,
                -1.0
            );

        assert(r.high == 0.0);
        assert(r.low == 0.0);
    }


    /*
     * estimate() is deliberately a rounded approximation.
     *
     * Components here follow expansion order: least significant first.
     */
    {
        ExpansionBuffer!4 expansion;

        expansion.append(1.0);
        expansion.append(2.0);
        expansion.append(4.0);

        assert(estimate(expansion) == 7.0);
    }


    /*
     * TwoSum: small addend lost from the rounded main result is
     * recovered exactly in the tail.
     */
    {
        const auto r =
            twoSum(
                10_000_000_000_000_000.0,
                1.0
            );

        assert(
            r.high ==
            10_000_000_000_000_000.0
        );

        assert(r.low == 1.0);
    }


    /*
     * Half an ulp at 1.0 rounds to even in the main component and is
     * retained exactly as the tail.
     */
    {
        enum double halfUlpAtOne =
            0x1p-53;

        const auto r =
            twoSum(
                1.0,
                halfUlpAtOne
            );

        assert(r.high == 1.0);
        assert(r.low == halfUlpAtOne);
    }


    /*
     * TwoDiff recovers a lost low-order unit.
     */
    {
        const auto r =
            twoDiff(
                10_000_000_000_000_000.0,
                1.0
            );

        assert(
            r.high ==
            10_000_000_000_000_000.0
        );

        assert(r.low == -1.0);
    }


    /*
     * Exact operations naturally produce a zero tail.
     */
    {
        const auto sum =
            twoSum(2.0, 4.0);

        assert(sum.high == 6.0);
        assert(sum.low == 0.0);

        const auto diff =
            twoDiff(7.0, 3.0);

        assert(diff.high == 4.0);
        assert(diff.low == 0.0);
    }


    /*
     * Split reconstructs representative binary64 values exactly after
     * ordinary binary64 addition.
     */
    {
        enum double[] values = [
            1.0,
            -1.0,
            0x1.23456789abcdep+100,
            -0x1.abcdef0123456p-100,
            0x1.0000000000001p+0
        ];

        static foreach (value; values)
        {{
            const auto parts =
                split(value);

            assert(
                roundedAdd(
                    parts.high,
                    parts.low
                ) == value
            );
        }}
    }


    /*
     * Known exact TwoProduct case.
     *
     * Let
     *
     *     a = 1 + 2^-52
     *
     * then
     *
     *     a² = 1 + 2^-51 + 2^-104
     *
     * The first two terms form the rounded binary64 high component;
     * 2^-104 remains as the exact tail.
     */
    {
        enum double a =
            0x1.0000000000001p+0;

        const auto r =
            twoProduct(a, a);

        assert(
            r.high ==
            0x1.0000000000002p+0
        );

        assert(
            r.low ==
            0x1p-104
        );
    }


    /*
     * Exact products have zero tails.
     */
    {
        const auto r =
            twoProduct(3.0, 4.0);

        assert(r.high == 12.0);
        assert(r.low == 0.0);
    }


    /*
     * Sign handling.
     */
    {
        const auto r =
            twoProduct(
                -0x1.0000000000001p+0,
                 0x1.0000000000001p+0
            );

        assert(
            r.high ==
            -0x1.0000000000002p+0
        );

        assert(
            r.low ==
            -0x1p-104
        );
    }


    /*
     * The helper types remain simple value types.
     */
    static assert(
        TwoComponent.sizeof ==
        2 * double.sizeof
    );

    static assert(
        SplitComponent.sizeof ==
        2 * double.sizeof
    );
}
