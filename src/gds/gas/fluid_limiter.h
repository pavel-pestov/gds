#ifndef _FLUID_CORRECTOR_H_
#define _FLUID_CORRECTOR_H_

#include "types.h"

namespace gds {

template<typename T>
Vector4<T> fluid(const Vector4<T>& x, T tdh)
{
    return Vector4<T>(x[0] * x[2] *tdh,
            (x[1] / (x[3] - 1.0) + x[0] * x[2] * x[2] * 0.5 + x[1]) * x[2] * tdh,
            (x[0] * x[2] * x[2] + x[1]) * tdh,
            x[3]);
}

template<typename T>
Vector4<T> fluid_limiter(const Vector4<T>& l, const Vector4<T>& m, const Vector4<T>& r, const T tdh)
{
    if (m[0] <= 0.0f)
        return Vector4<T>();
    Vector4<T> f = fluid(m, tdh);
    T cv = 1.0;
    T e;
    if (m[2] > 0.0)
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
    return Vector4<T>(m[0], m[1], m[2] * cv, m[3]);
}

}

#endif // _FLUID_CORRECTOR_H_
