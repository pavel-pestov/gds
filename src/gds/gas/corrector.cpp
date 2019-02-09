#include "corrector.h"

namespace gds {

template<typename T>
T order3(const T dl, const T dr, const T shift)
{
    T c1 = min(fabs(dl), fabs(dr));
    T c2 = c1 + fabs(dl - dr);
    T sh = 1 - (((dl < 0.0) ^ (dr < 0.0)) << 1);
    if (c1 * (3.0 + sh) < c2) {
        if (c1 * (5.0 + sh + sh) < c2) {
            sh = (2.0 + sh) * c1 / (6.0 * (c2 - c1)) ;
            c2 += c1;
            c1 /= c2 + c2;
            c2 = sh;
        } else {
            c1 = (1.25 - sh) * ((5.0 + 6.0 * sh) * c1 + c2) / (9.0 * (c1 + c2));
            c2 = 1.0 / 12.0;
        }
        if ((dl < 0.0) ^ (dr < 0.0)) {
            sh = c1;
            c1 = c2;
            c2 = sh;
        }
    } else {
        c1 = c2 = 1.0 / 12.0;
    }
    sh = 3.0 * (sign(shift) - shift);
    return c1 * (dr + dl) * sh + c2 * (dr - dl) * (1.0 - shift * (sh + shift));
}

template<typename T>
vector4<T> decomposition(const T dro, T dp, T dv, const T ro, const T bc, const T bc2)
{
    dp *= 0.5 * bc2;
    dv *= 0.5 * ro * bc;
    return vector4<T>(dp - dv, dro - dp - dp, dp + dv, 0);
}

template<typename T>
vector4<T> correction(vector4<T> l, vector4<T> x, vector4<T> r, T tdh)
{
    T c2 = x[3] * x[1] / x[0];
    T c = sqrt(c2);
    T bc = 1.0 / c;
    T bc2 = bc * bc;
    l = decomposition(x[0] - l[0], x[1] - l[1], x[2] - l[2], x[0], bc, bc2);
    r = decomposition(r[0] - x[0], r[1] - x[1], r[2] - x[2], x[0], bc, bc2);
    bc = tdh * x[2];
    tdh *= c;
    bc2 = x[0];
    x[0] = order3(l[0], r[0], bc - tdh);
    x[1] = order3(l[1], r[1], bc);
    x[2] = order3(l[2], r[2], bc + tdh);
    return vector4<T>((x[0] + x[2]) + x[1], c2 * (x[0] + x[2]), c * (x[2] - x[0]) / bc2, bc);
}

template<typename T>
vector4<T> corrector(const vector4<T>& l, const vector4<T>& x, const vector4<T>& r, const T tdh)
{
    float z = zero3(min(x[0], min(l[0], r[0])) * 1e3 / (max(x[0], max(l[0], r[0])) + 1e-20));
    if (z <= 0.0)
        return x;
    vector4<T> d = correction(l, x, r, tdh);
    if (!isfinite(d[0] + d[1] + d[2]))
        return x;
    d[0] *= z;
    d[1] *= z;
    d[2] *= z;
    if (d[0] <= -x[0] || d[1] <= -x[1])
        return x;
    return vector4<T>(x[0] + d[0], x[1] + d[1], x[2] + d[2], x[3]);
}

template<typename T>
T corrector(const T l, const T x, const T r, const T v, const T tdh)
{
    return x + high_order(x - l, r - x, v * tdh);
}

}
