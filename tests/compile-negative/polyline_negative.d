module polyline_negative;

import geo3;

@safe Polyline3View!double escape()
{
    Point3!double[2] local;

    return Polyline3View!double(local[]);
}
