#include "fluid_limiter.h"

namespace gds {

namespace math {

double4 fluid(const double4& x, double tdh)
{
    return double4(x[0] * x[2] *tdh,
            (x[1] / (x[3] - 1.0) + x[0] * x[2] * x[2] * 0.5 + x[1]) * x[2] * tdh,
            (x[0] * x[2] * x[2] + x[1]) * tdh,
            x[3]);
}

double E(const double4& x)
{
    return x[1] / (x[3] - 1.0) + x[0] * x[2] * x[2] * 0.5;
}

double4 fluid_corrector(const double4& l, const double4& x, const double4& r, const float tdh)
{
    if (x[0] <= 0.0f)
        return double4();
    double4 f = fluid(x, tdh);
    double cv = 1.0;
    double e;
    if (x[2] > 0.0)
    {
        e = E(l);
        if (l[0] < f[0])
            cv = l[0] / f[0];
        if (e < f[1] * cv)
            cv *= e / f[1];
    }
    else
    {
        e = E(r);
        if (r[0] < -f[0])
            cv = -r[0] / f[0];
        if (e < -f[1] * cv)
            cv *= -e / f[1];
    }
    return double4(x[0], x[1], x[2] * cv, x[3]);
}

FluidLimiter::FluidLimiter()
{

}

double4 FluidLimiter::operator()(const double4& l, const double4& x, const double4& r, const float tdh)
{
    return fluid_corrector(l, x, r, tdh);
}

}

}
