#ifndef _EULER_CORRECTOR_H_
#define _EULER_CORRECTOR_H_

#include "types.h"

namespace gds {

template<typename T>
T order3(T hl, T hm, T hr, T dl, T dr, T shift)
{
    T du = (dl + dr) * hm * 2.0 / (hl + hm + hm + hr);
    T d2u = ((hm + hl) * dr - (hm + hr) * dl) * hm * hm * 8.0 / ((hm + hl) * (hm + hr) * (hl + hm + hm + hr));
    shift /= hm;
    hl = fabs(du);
    hr = fabs(d2u);
    hm = fabs(dl) < fabs(dr);
    hm = fabs(hm * dl + (1.0 - hm) * dr);
    dr = !((dl < 0.0) ^ (dr <= 0.0));
    dl = hm * dr;
    dl = hr + 12.0 * dl < 6.0 * hl ? (dl + hr / 12.0) / hl : 0.5;
    if (((0.5 < dr && 4.0 * dl * hl < hr) || (dr < 0.5 && (hm - dl * hl) * 12.0 < hr)) && 0.0 < hl) {
        if (0.5 < dr) {
            dl = 1.5 * hm / hl;
            dr = dl * hl / (3.0 * hr);
        } else {
            dr = (hm - dl * hl + hr / 12.0) / (hr + hr);
            dl = (hm - dr * hr) / hl;
        }
    } else {
        dr = 1.0 / 12.0;
    }
    hm = 0.0 < shift;
    hm = hm + hm - 1.0;
    hl = hm - shift;
    hr = (hl + hl + hm) * shift - 1.0;
    return dl * du * hl - dr * d2u * hr;
}

template<typename T>
SimpleGas<T> decomposition(const T dro, T dp, T dv, const T ro, const T bc, const T bc2)
{
    dp *= 0.5 * bc2;
    dv *= 0.5 * ro * bc;
    return SimpleGas<T>(dp - dv, dro - dp - dp, dp + dv, 0);
}

template<typename T>
SimpleGas<T> correction(SimpleCell<T> l, SimpleCell<T> m, SimpleCell<T> r, T t)
{
    T c2 = m[3] * m[1] / m[0];
    T c = sqrt(c2);
    T bc = 1.0 / c;
    T bc2 = bc * bc;
    l.gas = decomposition(m[0] - l[0], m[1] - l[1], m[2] - l[2], m[0], bc, bc2);
    r.gas = decomposition(r[0] - m[0], r[1] - m[1], r[2] - m[2], m[0], bc, bc2);
    t / m.h;
    bc = t * m[2];
    t *= c;
    bc2 = m[0];
    m[0] = order3(l.h, m.h, r.h, l[0], r[0], bc - t);
    m[1] = order3(l.h, m.h, r.h, l[1], r[1], bc);
    m[2] = order3(l.h, m.h, r.h, l[2], r[2], bc + t);
    return SimpleGas<T>((m[0] + m[2]) + m[1], c2 * (m[0] + m[2]), c * (m[2] - m[0]) / bc2, bc);
}

template<typename T>
SimpleGas<T> corrector(const SimpleCell<T>& l, const SimpleCell<T>& m, const SimpleCell<T>& r, const T t)
{
    float z = zero3(min(m[0], min(l[0], r[0])) * 1e3 / (max(m[0], max(l[0], r[0])) + 1e-20));
    if (z <= 0.0)
        return m.gas;
    SimpleGas<T> d = correction(l, m, r, t);
    if (!isfinite(d[0] + d[1] + d[2]))
        return m.gas;
    d[0] *= z;
    d[1] *= z;
    d[2] *= z;
    if (d[0] <= -m[0] || d[1] <= -m[1])
        return m.gas;
    return SimpleGas<T>(m[0] + d[0], m[1] + d[1], m[2] + d[2], m[3]);
}

template<typename T>
T corrector(const T l, const T m, const T r, const T v, const T hl, const T hm, const T hr, const T t)
{
    return m + order3(hl, hm, hr, m - l, r - m, v * t);
}

}

#endif // _EULER_CORRECTOR_H_
