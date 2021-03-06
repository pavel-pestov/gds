#ifndef _REIMANN_SOLVER_H_
#define _REIMANN_SOLVER_H_

#include "types.h"

namespace gds {

template<typename T>
Vector4<T> vacum_solver(const Vector4<T>& g, T d)
{
    T c, cc, g1 =  2.0 / (g[3] - 1.0);
    if (g[0] <= 0.0)
        return Vector4<T>();
    if (g[1] <= 0.0)
        if ((d > 0) ? g[2] > 0 : g[2] < 0)
            return g;
        else
            return Vector4<T>();
    c = sqrt(g[3] * g[1] / g[0]);
    if (d > 0 ? g[2] <= -c * g1 : g[2] >= c * g1)
        return Vector4<T>();
    if (d > 0 ? g[2] >= c : g[2] <= -c)
        return g;
    cc = (c * 2.0 + d * (g[3] - 1.0) * g[2]) / (g[3] + 1.0);
    c = cc / c;
    return Vector4<T>(g[0] * pow(c, g1), g[1] * pow(c, g[3] * g1), d  * cc, g[3]);
}

template<typename T>
Vector4<T> safe_solver(const Vector4<T>& l, const Vector4<T>& r)
{
    T cl = sqrt(l[3] * l[1] / l[0]) + 1e-20;
    T cr = sqrt(r[3] * r[1] / r[0]) + 1e-20;
    T al = l[0] * (l[2] - min(l[2] - cl, r[2] - cr));
    T ar = r[0] * (max(l[2] + cl, r[2] + cr) - r[2]);
    T p = (l[1] * ar + r[1] * al - ar * al * (r[2] - l[2])) / (al + ar);
    T v = (r[2] * ar + l[2] * al - r[1] + l[1]) / (al + ar);
    if (!isfinite(p + v))
        return Vector4<T>(NAN);
    if (v >= 0) {
        if(l[2] - al / l[0] >= 0) {
            return l;
        } else {
            return Vector4<T>(1.0 / (1.0 / l[0] + (v - l[2]) / al), p, v, l[3]);
        }
    } else {
        if(r[2] + ar / r[0] <= 0) {
            return r;
        } else {
            return Vector4<T>(1.0 / (1.0 / r[0] - (v - r[2]) / ar), p, v, r[3]);
        }
    }
}

template<typename T>
Vector4<T> fast_solver(Vector4<T> l, Vector4<T> r)
{
    T cl = sqrt(l[3] * l[1] / l[0]);
    T cr = sqrt(r[3] * r[1] / r[0]);
    T al = l[0] * (l[2] - min(l[2] - cl, r[2] - cr));
    T ar = r[0] * (max(l[2] + cl, r[2] + cr) - r[2]);
    T p = (l[1] * ar + r[1] * al - ar * al * (r[2] - l[2])) / (al + ar);
    T v = (l[3] + 1.0) / (4.0 * l[3] * l[1]);
    T u = (r[3] + 1.0) / (4.0 * r[3] * r[1]);
    T ro = al;
    T n = ar;
    al = ro * (1.0 + (p - l[1]) * v);
    ar = n * (1.0 + (p - r[1]) * u);
    p = (l[1] * ar + r[1] * al - ar * al * (r[2] - l[2])) / (al + ar);
    al = (p - l[1]) * v;
    ar = (p - r[1]) * u;
    al = ro * (1.0 + al * (1.0 - 0.5 * al * (1.0 - al)));
    ar = n * (1.0 + ar * (1.0 - 0.5 * ar * (1.0 - ar)));
    v  = (al * l[2] + ar * r[2] - r[1] + l[1]) / (al + ar);
    if (!isfinite(p + v))
        return Vector4<T>(NAN);
    n = 1.0;
    if (v >= 0.0) {
        v = -v;
        r = l;
        r[2] = -r[2];
        cr = cl;
        n = -1.0;
    }
    al = 1.0 / (r[3] + 1.0);
    ar = 1.0 / r[1];
    l[1] = p * ar;
    cl = (r[3] - 1.0);
    if (p > r[1]) {
        l[2] = cl * al;
        l[2] = r[0] * (l[1] + l[2]) / (l[2] * l[1] + 1.0);
        if (v + (r[2] - v) * r[0] / (r[0] - l[2]) <= 0.0) {
            ro = r[0];
            p = r[1];
            v = r[2];
        } else {
            ro = l[2];
        }
    } else {
        l[0] = 1.0 / r[3];
        u = cl * 0.5;
        if (r[2] + cr <= 0.0) {
            ro = r[0];
            p = r[1];
            v = r[2];
        } else if (0.0 <= v + cr * pow(l[1], u * l[0])) {
            ro = r[0] * pow(l[1], l[0]);
        } else {
            l[2] = 2.0 * (cr - u * r[2]) * al;
            l[0] = l[2] / cr;
            cl = 2.0 / cl;
            ro = r[0] * pow(l[0], cl);
            p = r[1] * pow(l[0], cl * r[3]);
            v = -l[2];
        }
    }
    return Vector4<T>(ro, p, v * n, r[3]);
}

template<typename T>
Vector4<T> full_solver(Vector4<T> l, Vector4<T> r)
{
    T ro, p, v, f, d, dp;
    T cl = sqrt(l[3] * l[1] / l[0]);
    T cr = sqrt(r[3] * r[1] / r[0]);
    T al = l[1] > r[1] ? r[1] : l[1];
    T ar = l[1] > r[1] ? l[1] : r[1];
    p = 0.5 * (l[1] + r[1]) + 0.125 * (l[2] - r[2]) * (l[0] + r[0]) * (cl + cr);

    f = l[3] - 1.0;
    d = r[3] - 1.0;
    if (p < al) {
        ro = 1.0 / cr;
        v = 1.0 / cl;
        al = pow(l[1] / r[1], (f / l[3] + d / r[3]) * 0.25);
        ar = (al * l[2] * v + r[2] / cr + (1.0 / f + 1.0 / d) * (al - 1.0)) / (al * v + ro);
        al = 1.0 + f * 0.5 * (l[2] - ar) * v;
        ar = 1.0 + d * 0.5 * (ar - r[2]) * ro;
        p = 0.5 * (l[1] * pow(al, 2.0 * l[3] / f) + r[1] * pow(ar, 2.0 * r[3] / d));
    } else if (al < ar * 0.5 || p > ar) {
        al = sqrt(2.0 / (l[0] * (f * l[1] + (l[3] + 1.0) * p)));
        ar = sqrt(2.0 / (r[0] * (d * r[1] + (r[3] + 1.0) * p)));
        p = (al * l[1] + ar * r[1] - (r[2] - l[2])) / (al + ar);
    }

    al = (l[3] - 1.0) * l[1] / (l[3] + 1.0);
    ar = (r[3] - 1.0) * r[1] / (r[3] + 1.0);

    do {
        if (p <= l[1]) {
            v = p / l[1];
            f = 2.0 * cl * (pow(v, (l[3] - 1.0) / (2.0 * l[3])) - 1.0) / (l[3] - 1.0);
            d = pow(v, -(l[3] + 1.0) / (2.0 * l[3])) / (l[0] * cl);
        } else {
            v = sqrt(2.0 / ((l[3] + 1.0) * (al + p) * l[0]));
            f = (p - l[1]) * v;
            d = (1.0 - 0.5 * (p - l[1]) / (al + p)) * v;
        }
        if (p <= r[1]) {
            v = p / r[1];
            ro = 2.0 * cr * (pow(v, (r[3] - 1.0) / (2.0 * r[3])) - 1.0) / (r[3] - 1.0);
            d += pow(v, -(r[3] + 1.0) / (2.0 * r[3])) / (r[0] * cr);
        } else {
            v = sqrt(2.0 / ((r[3] + 1.0) * (ar + p) * r[0]));
            ro = (p - r[1]) * v;
            d += (1.0 - 0.5 * (p - r[1]) / (ar + p)) * v;
        }
        v = ro - f;
        f += ro;
        dp = (f + r[2] - l[2]) / d;
        p = fabs(p - dp);
    } while(fabs(dp) / p > 1e-3);

    v = 0.5 * (l[2] + r[2] + v);

    if (!isfinite(p+v))
        return Vector4<T>(NAN);

    dp = 1.0;
    if (v >= 0.0) {
        v = -v;
        r = l;
        r[2] = -r[2];
        cr = cl;
        dp = -1.0;
    }

    al = 1.0 / (r[3] + 1.0);
    ar = 1.0 / r[1];
    l[1] = p * ar;
    d = (r[3] - 1.0);
    if (p > r[1]) {
        l[2] = d * al;
        l[2] = r[0] * (l[1] + l[2]) / (l[2] * l[1] + 1.0);
        if (v + (r[2] - v) * r[0] / (r[0] - l[2]) <= 0.0) {
            ro = r[0];
            p = r[1];
            v = r[2];
        } else {
            ro = l[2];
        }
    } else {
        f = 1.0 / r[3];
        if (r[2] + cr <= 0.0) {
            ro = r[0];
            p = r[1];
            v = r[2];
        } else if (0.0 <= v + cr * pow(l[1], d * 0.5 * f)) {
            ro = r[0] * pow(l[1], f);
        } else {
            l[2] = 2.0 * (cr - d * 0.5 * r[2]) * al;
            f = l[2] / cr;
            d = 2.0 / d;
            ro = r[0] * pow(f, d);
            p = r[1] * pow(f, d * r[3]);
            v = -l[2];
        }
    }
    return Vector4<T>(ro, p, v * dp, r[3]);
}

template<typename T>
Vector4<T> rieman_solver(const Vector4<T>& l, const Vector4<T>& r)
{
    Vector4<T> fluid(NAN);
    T cl = sqrt(l[1] * l[3] / l[0]);
    T cr = sqrt(r[1] * r[3] / r[0]);
    if (r[0] <= 1e-20 || cr * (2.0 / (r[3] - 1.0)) < r[2]) {
        fluid = vacum_solver(l, 1.0);
    } else if (l[0] <= 1e-20 || -cl * (2.0 / (l[3] - 1.0)) > l[2]) {
        fluid = vacum_solver(r, -1.0);
    } else if (r[1] * l[1] > 0.0) {
        cl = (l[2] - r[2]) / (cl + cr);
        cr = fabs(l[1] - r[1]) / (l[1] + r[1]);
        if(0.04 * (cr * cr * cr * 2.0 + cl * cl) < 0.001) {
            fluid = fast_solver(l, r);
        } else {
            fluid = full_solver(l, r);
        }
    }
    if (!isfinite(fluid[0] + fluid[1] + fluid[2] + fluid[3])) {
        fluid = safe_solver(l, r);
    }
    return fluid;
}

}

#endif // _REIMANN_SOLVER_H_
