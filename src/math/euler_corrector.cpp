#include "euler_corrector.h"

namespace gds {

namespace math {

double correction(const double dl, const double dr, const double shift)
{
    double sh;
    double c1 = fabs(dr);
    double c2 = fabs(dl);
    if (c2 < c1)
    {
        sh = c1;
        c1 = c2;
        c2 = sh;
    }
    sh = dl * dr;
    if (sh > 0.0 && c1 * 4.0 < c2) {
        sh = c1 + c2;
        if (c1 * 8.0 < sh) {
            c2 = 0.5 * c1 / (c2 - c1);
            c1 = 1.5 * c1 / sh;
        } else {
            c1 = (c1 * 10.0 + sh) / (sh * 12.0);
            c2 = 1.0 / 12.0;
        }
    } else if (sh <= 0.0 && c1 * 3.5 < c2) {
        sh = c2 - c1;
        if (c1 * 4.0 <= sh) {
            c1 = c1 / sh;
            c2 = 0.0;
        } else {
            c2 = (c1 - sh * 0.25) / (c1 + c2);
            c1 = 0.25;
        }
    } else {
        c1 = 0.25;
        c2 = 1.0 / 12.0;
    }
    sh = ((shift < 0.0) ? -1.0 : 1.0);
    return c1 * (dr + dl) * (sh - shift) + c2 * (dr - dl) * (shift * (shift * 2.0 - sh * 3.0) + 1.0);
}

double4 decomposition(const double dro, double dp, double dv, const double ro, const double bc, const double bc2)
{
    dp *= 0.5 * bc2;
    dv *= 0.5 * ro * bc;
    return double4(dp - dv, dro - dp - dp, dp + dv, 0);
}

double4 gas_correction(double4 l, double4 x, double4 r, double tdh)
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
    x[0] = correction(l[0], r[0], bc - tdh);
    x[1] = correction(l[1], r[1], bc);
    x[2] = correction(l[2], r[2], bc + tdh);
    return double4((x[0] + x[2]) + x[1], c2 * (x[0] + x[2]), c * (x[2] - x[0]) / bc2, bc);
}

double4 euler_correction(const double4& l, const double4& x, const double4& r, const double tdh)
{
    float z = zero3(min(x[0], min(l[0], r[0])) * 1e3 / (max(x[0], max(l[0], r[0])) + 1e-20));
    if (z <= 0.0)
        return x;
    double4 d = gas_correction(l, x, r, tdh);
    if (!isfinite(d[0] + d[1] + d[2]))
        return x;
    d[0] *= z;
    d[1] *= z;
    d[2] *= z;
    if (d[0] <= -x[0] || d[1] <= -x[1])
        return x;
    return double4(x[0] + d[0], x[1] + d[1], x[2] + d[2], x[3]);
}

EulerCorrector::EulerCorrector()
{

}

double4 EulerCorrector::operator()(const double4& l, const double4& x, const double4& r, const double tdh)
{
    return euler_correction(l, x, r, tdh);
}

double EulerCorrector::operator()(const double l, const double x, const double r, const double v, const double tdh)
{
    return x + correction(x - l, r - x, v * tdh);
}

}

}
