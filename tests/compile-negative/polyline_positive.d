module polyline_positive;

import geo3;

@safe void probe()
{
    Point3!double[2] local;

    auto view =
        Polyline3View!double(local[]);

    assert(view.length == 2);
}
