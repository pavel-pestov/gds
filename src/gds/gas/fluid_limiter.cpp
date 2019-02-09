#include "fluid_limiter.h"

namespace gds {

template<typename T>
vector4<T> fluid(const vector4<T>& x, T tdh)
{
    return vector4<T>(x[0] * x[2] *tdh,
            (x[1] / (x[3] - 1.0) + x[0] * x[2] * x[2] * 0.5 + x[1]) * x[2] * tdh,
            (x[0] * x[2] * x[2] + x[1]) * tdh,
            x[3]);
}

template<typename T>
vector4<T> fluid_limiter(const vector4<T>& l, const vector4<T>& x, const vector4<T>& r, const T tdh)
{
    if (x[0] <= 0.0f)
        return vector4<T>();
    vector4<T> f = fluid(x, tdh);
    T cv = 1.0;
    T e;
    if (x[2] > 0.0)
    {
        e = l.e();
        if (l[0] < f[0])
            cv = l[0] / f[0];
        if (e < f[1] * cv)
            cv *= e / f[1];
    }
    else
    {
        e =r.e();
        if (r[0] < -f[0])
            cv = -r[0] / f[0];
        if (e < -f[1] * cv)
            cv *= -e / f[1];
    }
    return vector4<T>(x[0], x[1], x[2] * cv, x[3]);
}

}
