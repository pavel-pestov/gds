#ifndef _EULER_SOLVER_H_
#define _EULER_SOLVER_H_

#include "types.h"

namespace gds {

namespace math {

class EulerSolver
{
public:
    EulerSolver();
    double4 operator()(const double4& l, const double4& x, const double4& r, const float tdh);
};

}

}

#endif // _EULER_SOLVER_H_
