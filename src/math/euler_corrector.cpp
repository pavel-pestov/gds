#include "euler_corrector.h"

namespace gds {

namespace math {

double EulerCorrector::third_order_plus(const double dl, const double dr, const double shift)
{
    double sh;
    double c1 = fabs(dr);
    double c2 = fabs(dl);
    if (c2 < c1) {
        sh = c1;
        c1 = c2;
        c2 = sh;
    }
    sh = sign(dl) * sign(dr);
    if (c1 * (3.0 + sh) < c2) {
        sh /= (c1 + c2);
        if (sh > 0) {
            if (c1 * 7.0 < c2) {
                sh *= 1.5 * c1;
                c2 = 0.5 * c1 / (c2 - c1);
                c1 = sh;
            } else {
                c1 = (1.0 / 12.0) * (c1 * 11.0 + c2) * sh;
                c2 = 1.0 / 12.0;
            }
        } else {
            if (c1 * 3.0 < c2) {
                sh *= -0.5 * c1;
                c1 = 0.5 * c1 / (c2 - c1);
                c2 = sh;
            } else {
                c2 = -0.25 * (c2 - c1) * sh;
                c1 = 0.25;
            }
        }
    } else {
        c1 = 0.25;
        c2 = 1.0 / 12.0;
    }
    sh = sign(shift);
    return c1 * (dr + dl) * (sh - shift) + c2 * (dr - dl) * (shift * (shift + shift - sh * 3.0) + 1.0);
}

double EulerCorrector::third_order(const double dl, const double dr, const double shift)
{
    double c = fabs(dl - dr);
    double sh = min(fabs(dl), fabs(dr));
    c = (sign(dl) == sign(dr) ? -c : c) + fabs(dl + dr) * 3.0;
    c = sh * 12.0 < c ? sh / c : 1.0 / 12.0;
    sh = sign(shift);
    return c * ((dr + dl) * (sh - shift) * 3.0 + (dr - dl) * (shift * (shift + shift - sh * 3.0) + 1.0));
}

double EulerCorrector::second_order_plus(const double dl, const double dr, const double shift)
{
    double c = min(fabs(dl), fabs(dr));
    c = c * 4.0 < fabs(dl + dr) ? c / fabs(dl + dr) : 0.25;
    return c * (sign(shift) - shift) * (dl + dr);
}

double EulerCorrector::second_order(const double dl, const double dr, const double shift)
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
    case Order::second:
        high_order = &second_order;
        break;
    case Order::second_plus:
        high_order = &second_order_plus;
        break;
    case Order::third:
        high_order = &third_order;
        break;
    case Order::third_plus:
        high_order = &third_order_plus;
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
