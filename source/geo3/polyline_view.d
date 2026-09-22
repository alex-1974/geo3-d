/**
 * Non-owning views of ordered three-dimensional point sequences.
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
module geo3.polyline_view;

import geo3.point : Point3;
import geo3.scalar : isGeoScalar;
import geo3.segment : Segment3;


/**
 * Non-owning read-only view of an ordered sequence of 3D points.
 *
 * Supported scalar types are `int`, `long`, `float`, `double`, and `real`.
 *
 * `Polyline3View.init` is an empty view.
 *
 * Polyline3View does not allocate or copy point data. The caller retains
 * ownership of the backing storage, which must remain valid for the
 * lifetime of the view.
 *
 * The view aliases its backing storage. Changes made to mutable backing
 * storage through its owner remain visible through an existing view.
 * Mutation is not exposed through Polyline3View itself.
 *
 * Empty and singleton polylines are valid.
 */
struct Polyline3View(T)
if (isGeoScalar!T)
{
private:
    const(Point3!T)[] _points;

public:
    /**
     * Constructs a view over contiguous point storage.
     *
     * No point data is copied.
     */
    this(return scope const(Point3!T)[] points)
        pure nothrow @safe @nogc
    {
        _points = points;
    }


    /// Number of stored points.
    @property size_t length() const
        pure nothrow @safe @nogc
    {
        return _points.length;
    }


    /// True when the view contains no points.
    @property bool empty() const
        pure nothrow @safe @nogc
    {
        return _points.length == 0;
    }


    /**
     * Number of consecutive line segments represented by the view.
     *
     * Empty and singleton polylines contain no segments.
     */
    @property size_t segmentCount() const
        pure nothrow @safe @nogc
    {
        return _points.length > 1
            ? _points.length - 1
            : 0;
    }


    /**
     * Returns the segment beginning at point index.
     *
     * Valid indices are:
     *
     *     0 .. segmentCount
     */
    Segment3!T segment(size_t index) const
        pure nothrow @safe @nogc
    {
        return Segment3!T(
            _points[index],
            _points[index + 1]
        );
    }


    /**
     * Returns one point by value.
     *
     * Mutation of the backing storage is not exposed through the view.
     */
    Point3!T opIndex(size_t index) const
        pure nothrow @safe @nogc
    {
        return Point3!T(
            _points[index].x,
            _points[index].y,
            _points[index].z
        );
    }
}


/// Example borrowing caller-owned point storage as a polyline.
@safe unittest
{
    import geo3;

    alias P = Point3!double;
    alias S = Segment3!double;
    alias V = Polyline3View!double;

    P[3] points = [
        P(0.0, 0.0, 0.0),
        P(2.0, 1.0, 3.0),
        P(5.0, 1.0, 7.0)
    ];

    auto polyline =
        V(points[]);

    assert(!polyline.empty);
    assert(polyline.length == 3);
    assert(polyline.segmentCount == 2);

    assert(
        polyline[1] ==
        P(2.0, 1.0, 3.0)
    );

    assert(
        polyline.segment(1) ==
        S(
            P(2.0, 1.0, 3.0),
            P(5.0, 1.0, 7.0)
        )
    );

    /*
     * The view aliases rather than copies caller-owned storage.
     */
    points[1] =
        P(3.0, 2.0, 4.0);

    assert(
        polyline[1] ==
        P(3.0, 2.0, 4.0)
    );
}


@safe unittest
{
    import std.meta : AliasSeq;


    /*
     * Polyline3View follows the Point3 scalar domain.
     */
    static foreach (T; AliasSeq!(int, long, float, double, real))
    {
        static assert(Polyline3View!T.init.length == 0);
        static assert(Polyline3View!T.init.empty);
        static assert(Polyline3View!T.init.segmentCount == 0);
    }


    /*
     * Unsupported scalar types remain unavailable.
     */
    static assert(!__traits(compiles, Polyline3View!byte));
    static assert(!__traits(compiles, Polyline3View!short));
    static assert(!__traits(compiles, Polyline3View!uint));
    static assert(!__traits(compiles, Polyline3View!ulong));


    alias P = Point3!double;
    alias S = Segment3!double;
    alias V = Polyline3View!double;


    /*
     * Empty storage produces an empty view.
     */
    P[] noPoints;

    auto emptyView =
        V(noPoints);

    assert(emptyView.empty);
    assert(emptyView.length == 0);
    assert(emptyView.segmentCount == 0);


    /*
     * A singleton is valid and contains no segment.
     */
    P[1] singletonPoints = [
        P(1.0, 2.0, 3.0)
    ];

    auto singleton =
        V(singletonPoints[]);

    assert(!singleton.empty);
    assert(singleton.length == 1);
    assert(singleton.segmentCount == 0);
    assert(singleton[0] == P(1.0, 2.0, 3.0));


    /*
     * Multi-point views preserve point order and all three coordinates.
     */
    P[3] points = [
        P(1.0, 2.0, 3.0),
        P(4.0, 5.0, 6.0),
        P(7.0, 8.0, 9.0)
    ];

    auto view =
        V(points[]);

    assert(view.length == 3);
    assert(view.segmentCount == 2);

    assert(view[0] == points[0]);
    assert(view[1] == points[1]);
    assert(view[2] == points[2]);

    assert(
        view.segment(0) ==
        S(
            P(1.0, 2.0, 3.0),
            P(4.0, 5.0, 6.0)
        )
    );

    assert(
        view.segment(1) ==
        S(
            P(4.0, 5.0, 6.0),
            P(7.0, 8.0, 9.0)
        )
    );


    /*
     * The view does not own or copy the backing storage.
     *
     * Mutation through the owner remains visible through an existing
     * read-only view, including changes to the Z coordinate.
     */
    points[1] =
        P(40.0, 50.0, 60.0);

    assert(
        view[1] ==
        P(40.0, 50.0, 60.0)
    );


    /*
     * Mutation is not available through Polyline3View.
     */
    static assert(
        !__traits(
            compiles,
            {
                P[1] backing;

                auto readOnly =
                    V(backing[]);

                readOnly[0] =
                    P(7.0, 8.0, 9.0);
            }
        )
    );


    /*
     * Immutable backing storage is directly viewable.
     */
    immutable P[2] immutablePoints = [
        P(-1.0, 2.0, 3.0),
        P(4.0, -5.0, 6.0)
    ];

    auto immutableView =
        V(immutablePoints[]);

    assert(immutableView.length == 2);
    assert(immutableView.segmentCount == 1);
    assert(immutableView[0] == P(-1.0, 2.0, 3.0));
    assert(immutableView[1] == P(4.0, -5.0, 6.0));


    /*
     * DIP1000 must reject a view escaping stack-owned backing storage.
     */
    static assert(
        !__traits(
            compiles,
            {
                @safe V invalidEscape()
                {
                    P[2] local;

                    return V(local[]);
                }
            }
        )
    );
}
