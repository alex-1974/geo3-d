/**
 * Non-owning views of implicitly closed three-dimensional linear rings.
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
module geo3.linear_ring_view;

import geo3.point : Point3;
import geo3.scalar : isGeoScalar;
import geo3.segment : Segment3;


/**
 * Non-owning read-only view of an ordered cyclic sequence of 3D points.
 *
 * Supported scalar types are `int`, `long`, `float`, `double`, and `real`.
 *
 * `LinearRing3View.init` is an empty ring.
 *
 * LinearRing3View does not allocate or copy point data. The caller retains
 * ownership of the backing storage, which must remain valid for the
 * lifetime of the view.
 *
 * The view aliases its backing storage. Changes made to mutable backing
 * storage through its owner remain visible through an existing view.
 * Mutation is not exposed through LinearRing3View itself.
 *
 * Closure is implicit: for a non-empty ring the final stored vertex is
 * connected back to the first stored vertex.
 *
 * No normalization is performed. An explicitly repeated first vertex
 * remains an ordinary stored vertex and participates in the implicit
 * closing traversal.
 *
 * No planarity requirement is imposed by this representation. A non-planar
 * cyclic point sequence is representable; geometric validation is a
 * separate operation.
 *
 * Empty and degenerate rings are valid representations.
 */
struct LinearRing3View(T)
if (isGeoScalar!T)
{
private:
    const(Point3!T)[] _points;

public:
    /**
     * Constructs a view over contiguous point storage.
     *
     * No point data is copied and no normalization is performed.
     */
    this(return scope const(Point3!T)[] points)
        pure nothrow @safe @nogc
    {
        _points = points;
    }


    /// Number of stored vertices.
    @property size_t length() const
        pure nothrow @safe @nogc
    {
        return _points.length;
    }


    /// True when the ring contains no stored vertices.
    @property bool empty() const
        pure nothrow @safe @nogc
    {
        return _points.length == 0;
    }


    /**
     * Number of segments in the cyclic traversal.
     *
     * An empty ring contains no segments. Every non-empty ring contains
     * one segment per stored vertex, including the implicit closing
     * segment.
     */
    @property size_t segmentCount() const
        pure nothrow @safe @nogc
    {
        return _points.length;
    }


    /**
     * Returns one vertex by value.
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


    /**
     * Returns one segment in stored traversal order.
     *
     * The final segment closes the ring by connecting the final stored
     * vertex back to the first stored vertex.
     *
     * Valid indices are:
     *
     *     0 .. segmentCount
     */
    Segment3!T segment(size_t index) const
        pure nothrow @safe @nogc
    {
        /*
         * Access first so an invalid index retains the normal D bounds
         * semantics before calculating the cyclic successor.
         */
        const current = _points[index];

        const size_t next =
            index + 1 == _points.length
                ? 0
                : index + 1;

        return Segment3!T(
            Point3!T(
                current.x,
                current.y,
                current.z
            ),
            Point3!T(
                _points[next].x,
                _points[next].y,
                _points[next].z
            )
        );
    }
}


/// Example using the implicit closing edge of a linear ring.
@safe unittest
{
    import geo3;

    alias P = Point3!double;
    alias S = Segment3!double;
    alias R = LinearRing3View!double;

    P[3] points = [
        P(0.0, 0.0, 1.0),
        P(4.0, 0.0, 2.0),
        P(0.0, 3.0, 3.0)
    ];

    auto ring =
        R(points[]);

    assert(!ring.empty);
    assert(ring.length == 3);
    assert(ring.segmentCount == 3);

    assert(
        ring[1] ==
        P(4.0, 0.0, 2.0)
    );

    /*
     * Ring closure is implicit: the final segment returns to vertex 0.
     */
    assert(
        ring.segment(2) ==
        S(
            P(0.0, 3.0, 3.0),
            P(0.0, 0.0, 1.0)
        )
    );
}


@safe unittest
{
    import std.meta : AliasSeq;


    /*
     * LinearRing3View follows the Point3 scalar domain.
     */
    static foreach (T; AliasSeq!(int, long, float, double, real))
    {
        static assert(LinearRing3View!T.init.length == 0);
        static assert(LinearRing3View!T.init.empty);
        static assert(LinearRing3View!T.init.segmentCount == 0);
    }


    /*
     * Unsupported scalar types remain unavailable.
     */
    static assert(!__traits(compiles, LinearRing3View!byte));
    static assert(!__traits(compiles, LinearRing3View!short));
    static assert(!__traits(compiles, LinearRing3View!uint));
    static assert(!__traits(compiles, LinearRing3View!ulong));


    alias P = Point3!double;
    alias S = Segment3!double;
    alias R = LinearRing3View!double;


    /*
     * Empty storage produces an empty ring with no segments.
     */
    P[] noPoints;

    auto emptyRing =
        R(noPoints);

    assert(emptyRing.empty);
    assert(emptyRing.length == 0);
    assert(emptyRing.segmentCount == 0);


    /*
     * A singleton is a valid degenerate ring.
     *
     * Its only segment closes the sole vertex onto itself.
     */
    P[1] singletonPoints = [
        P(1.0, 2.0, 3.0)
    ];

    auto singleton =
        R(singletonPoints[]);

    assert(!singleton.empty);
    assert(singleton.length == 1);
    assert(singleton.segmentCount == 1);
    assert(singleton[0] == P(1.0, 2.0, 3.0));

    assert(
        singleton.segment(0) ==
        S(
            P(1.0, 2.0, 3.0),
            P(1.0, 2.0, 3.0)
        )
    );


    /*
     * Two vertices produce two opposite directed segments.
     */
    P[2] twoPoints = [
        P(1.0, 2.0, 3.0),
        P(4.0, 5.0, 6.0)
    ];

    auto twoVertexRing =
        R(twoPoints[]);

    assert(twoVertexRing.length == 2);
    assert(twoVertexRing.segmentCount == 2);

    assert(
        twoVertexRing.segment(0) ==
        S(
            P(1.0, 2.0, 3.0),
            P(4.0, 5.0, 6.0)
        )
    );

    assert(
        twoVertexRing.segment(1) ==
        S(
            P(4.0, 5.0, 6.0),
            P(1.0, 2.0, 3.0)
        )
    );


    /*
     * Ordinary rings preserve stored vertex order and close implicitly.
     */
    P[3] points = [
        P(0.0, 0.0, 1.0),
        P(4.0, 0.0, 2.0),
        P(0.0, 3.0, 3.0)
    ];

    auto ring =
        R(points[]);

    assert(ring.length == 3);
    assert(ring.segmentCount == 3);

    assert(ring[0] == points[0]);
    assert(ring[1] == points[1]);
    assert(ring[2] == points[2]);

    assert(
        ring.segment(0) ==
        S(points[0], points[1])
    );

    assert(
        ring.segment(1) ==
        S(points[1], points[2])
    );

    assert(
        ring.segment(2) ==
        S(points[2], points[0])
    );


    /*
     * Non-planar point sequences remain valid representations.
     *
     * Planarity is not a LinearRing3View invariant.
     */
    P[4] nonPlanarPoints = [
        P(0.0, 0.0, 0.0),
        P(4.0, 0.0, 0.0),
        P(4.0, 4.0, 1.0),
        P(0.0, 4.0, 0.0)
    ];

    auto nonPlanar =
        R(nonPlanarPoints[]);

    assert(nonPlanar.length == 4);
    assert(nonPlanar.segmentCount == 4);
    assert(nonPlanar[2] == P(4.0, 4.0, 1.0));

    assert(
        nonPlanar.segment(3) ==
        S(
            P(0.0, 4.0, 0.0),
            P(0.0, 0.0, 0.0)
        )
    );


    /*
     * A repeated final vertex is retained as supplied.
     *
     * It is not interpreted or removed as an explicit closure marker.
     */
    P[4] explicitlyClosed = [
        P(0.0, 0.0, 1.0),
        P(4.0, 0.0, 2.0),
        P(0.0, 3.0, 3.0),
        P(0.0, 0.0, 1.0)
    ];

    auto repeated =
        R(explicitlyClosed[]);

    assert(repeated.length == 4);
    assert(repeated.segmentCount == 4);

    assert(
        repeated.segment(2) ==
        S(
            P(0.0, 3.0, 3.0),
            P(0.0, 0.0, 1.0)
        )
    );

    assert(
        repeated.segment(3) ==
        S(
            P(0.0, 0.0, 1.0),
            P(0.0, 0.0, 1.0)
        )
    );


    /*
     * The ring view aliases rather than copies mutable backing storage.
     *
     * Z-coordinate changes through the owner remain visible.
     */
    points[1] =
        P(8.0, 0.0, 20.0);

    assert(
        ring[1] ==
        P(8.0, 0.0, 20.0)
    );

    assert(
        ring.segment(0) ==
        S(
            P(0.0, 0.0, 1.0),
            P(8.0, 0.0, 20.0)
        )
    );


    /*
     * Mutation is not available through LinearRing3View.
     */
    static assert(
        !__traits(
            compiles,
            {
                P[1] backing;

                auto readOnly =
                    R(backing[]);

                readOnly[0] =
                    P(7.0, 8.0, 9.0);
            }
        )
    );


    /*
     * Immutable backing storage is directly viewable.
     */
    immutable P[3] immutablePoints = [
        P(-1.0, 0.0, 1.0),
        P(1.0, 0.0, 2.0),
        P(0.0, 2.0, 3.0)
    ];

    auto immutableRing =
        R(immutablePoints[]);

    assert(immutableRing.length == 3);
    assert(immutableRing.segmentCount == 3);
    assert(
        immutableRing[0] ==
        P(-1.0, 0.0, 1.0)
    );


    /*
     * DIP1000 must reject a view escaping stack-owned backing storage.
     */
    static assert(
        !__traits(
            compiles,
            {
                @safe R invalidEscape()
                {
                    P[3] local;

                    return R(local[]);
                }
            }
        )
    );
}
