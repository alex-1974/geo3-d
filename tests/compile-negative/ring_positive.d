module ring_positive;

import geo3;

@safe void probe()
{
    Point3!double[3] local;

    auto ring =
        LinearRing3View!double(local[]);

    assert(ring.length == 3);
}
