#include "euler_solver.h"

namespace gds {

namespace math {

double4 euler_solver(const double4& l, const double4& x, const double4& r, const float tdh)
{
    double4 n;
    n[0] = x[0] + (l[0] * l[2] - r[0] * r[2]) * tdh;
    if (n[0] <= 1e-20)
        return double4();
    n[2] = 1.0f / n[0];
    n[3] = n[2] * tdh;
    n[2] = x[0] * x[2] * n[2];
    n[2] += (l[0] * l[2] * l[2] - r[0] * r[2] * r[2]) * n[3];
    n[2] += (l[1] - r[1]) * n[3];
    n[1] = x[1] / (x[3] - 1.0);
    n[3] = x[0] * x[2] * x[2] - n[0] * n[2] * n[2];
    n[1] += (l[1] * l[2] * l[3] / (l[3] - 1.0) - r[1] * r[2] * r[3] / (r[3] - 1.0)) * tdh;
    n[3] += (l[0] * l[2] * l[2] * l[2] - r[0] * r[2] * r[2] * r[2]) * tdh;
    n[1] = (n[1] + n[3] * 0.5) * (x[3] - 1.0);
    if (n[1] < 0.0)
        n[1] = 0.0;
    n[3] = x[3];
    return n;
}

EulerSolver::EulerSolver()
{

}

double4 EulerSolver::operator()(const double4& l, const double4& x, const double4& r, const float tdh)
{
    return euler_solver(l, x, r, tdh);
}

}

}
