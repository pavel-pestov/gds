#ifndef _EULER_SOLVER_H_
#define _EULER_SOLVER_H_

#include "types.h"

namespace gds {

template<typename T>
Vector4<T> euler_limiter(Vector4<T> x)
{
    if (fabs(x[2]) > sqrt((fabs(x[1] * x[3]) + 1e-20) / (fabs(x[0]) + 1e-20)))
        x[2] *= 0.5;
    else if (x[0] < x[1])
        x[1] = (fabs(x[1]) + 1.0) * 0.5;
    else
        x[0] = (fabs(x[0]) + 1.0) * 0.5;
    return x;
}

template<typename T>
Vector4<T> euler_solver(const Vector4<T>& l, const Vector4<T>& m, const Vector4<T>& r, const T tdh)
{
    Vector4<T> n;
    n[0] = m[0] + (l[0] * l[2] - r[0] * r[2]) * tdh;
    if (n[0] <= 1e-20)
        return Vector4<T>();
    n[2] = 1.0f / n[0];
    n[3] = n[2] * tdh;
    n[2] = m[0] * m[2] * n[2];
    n[2] += (l[0] * l[2] * l[2] - r[0] * r[2] * r[2]) * n[3];
    n[2] += (l[1] - r[1]) * n[3];
    n[1] = m[1] / (m[3] - 1.0);
    n[3] = m[0] * m[2] * m[2] - n[0] * n[2] * n[2];
    n[1] += (l[1] * l[2] * l[3] / (l[3] - 1.0) - r[1] * r[2] * r[3] / (r[3] - 1.0)) * tdh;
    n[3] += (l[0] * l[2] * l[2] * l[2] - r[0] * r[2] * r[2] * r[2]) * tdh;
    n[1] = (n[1] + n[3] * 0.5) * (m[3] - 1.0);
    if (n[1] < 0.0)
        n[1] = 0.0;
    n[3] = m[3];
    if (!isfinite(n[0] + n[1] + n[2] + n[3]))
        n = euler_limiter(m);
    return n;
}

}

#endif // _EULER_SOLVER_H_
