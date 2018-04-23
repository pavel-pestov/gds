#include "euler_corrector.h"

namespace gds {

namespace math {

double EulerCorrector::third_order_full(const double dl, const double dr, const double shift)
{
    double c1 = min(fabs(dl), fabs(dr));
    double c2 = c1 + fabs(dl - dr);
    double sh = 1 - (((dl < 0.0) ^ (dr < 0.0)) << 1);
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
        if ((dl < 0.0) ^ (dr < 0.0))
        {
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

double EulerCorrector::third_order_fast(const double dl, const double dr, const double shift)
{
    double c = fabs(dl - dr);
    double sh = 1 - (((dl < 0.0) ^ (dr < 0.0)) << 1);
    c = 3.0 * fabs(dl + dr) - sh * c;
    sh = min(fabs(dl), fabs(dr));
    c = 12.0 * sh < c ? sh / c : 1.0 / 12.0;
    sh = 3.0 * (sign(shift) - shift);
    return c * ((dr + dl) * sh + (dr - dl) * (1.0 - shift * (sh + shift)));
}

double EulerCorrector::second_order_full(const double dl, const double dr, const double shift)
{
    double c = min(fabs(dl), fabs(dr));
    double sh = fabs(dl + dr);
    c = 4.0 * c < sh ? c / sh : 0.25;
    sh = sign(shift) - shift;
    return c * sh * (dl + dr);
}

double EulerCorrector::second_order_fast(const double dl, const double dr, const double shift)
{
    double c = dl * dr;
    return c <= 0.0 ? 0.0 : (sign(shift) - shift) * c / (dl + dr);
}

double EulerCorrector::first_order(const double dl, const double dr, const double shift)
{
    return 0.0;
}


double4 EulerCorrector::decomposition(const double dro, double dp, double dv, const double ro, const double bc, const double bc2)
{
    dp *= 0.5 * bc2;
    dv *= 0.5 * ro * bc;
    return double4(dp - dv, dro - dp - dp, dp + dv, 0);
}

double4 EulerCorrector::correction(double4 l, double4 x, double4 r, double tdh)
{
    double c2 = x[3] * x[1] / x[0];
    double c = sqrt(c2);
    double bc = 1.0 / c;
    double bc2 = bc * bc;
    l = decomposition(x[0] - l[0], x[1] - l[1], x[2] - l[2], x[0], bc, bc2);
    r = decomposition(r[0] - x[0], r[1] - x[1], r[2] - x[2], x[0], bc, bc2);
    bc = tdh * x[2];
    tdh *= c;
    bc2 = x[0];
    x[0] = high_order(l[0], r[0], bc - tdh);
    x[1] = high_order(l[1], r[1], bc);
    x[2] = high_order(l[2], r[2], bc + tdh);
    return double4((x[0] + x[2]) + x[1], c2 * (x[0] + x[2]), c * (x[2] - x[0]) / bc2, bc);
}

double4 EulerCorrector::euler_correction(const double4& l, const double4& x, const double4& r, const double tdh)
{
    float z = zero3(min(x[0], min(l[0], r[0])) * 1e3 / (max(x[0], max(l[0], r[0])) + 1e-20));
    if (z <= 0.0)
        return x;
    double4 d = correction(l, x, r, tdh);
    if (!isfinite(d[0] + d[1] + d[2]))
        return x;
    d[0] *= z;
    d[1] *= z;
    d[2] *= z;
    if (d[0] <= -x[0] || d[1] <= -x[1])
        return x;
    return double4(x[0] + d[0], x[1] + d[1], x[2] + d[2], x[3]);
}

EulerCorrector::EulerCorrector(const Order order)
{
    switch (order) {
    case Order::second_fast:
        high_order = &second_order_fast;
        break;
    case Order::second_full:
        high_order = &second_order_full;
        break;
    case Order::third_fast:
        high_order = &third_order_fast;
        break;
    case Order::third_full:
        high_order = &third_order_full;
        break;
    default:
        high_order = &first_order;
        break;
    }
}

double4 EulerCorrector::operator()(const double4& l, const double4& x, const double4& r, const double tdh)
{
    return euler_correction(l, x, r, tdh);
}

double EulerCorrector::operator()(const double l, const double x, const double r, const double v, const double tdh)
{
    return x + high_order(x - l, r - x, v * tdh);
}

}

}
