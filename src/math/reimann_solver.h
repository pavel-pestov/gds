#ifndef _REIMANN_SOLVER_H_
#define _REIMANN_SOLVER_H_

#include "types.h"

namespace gds {

namespace math {

class ReimannSolver
{
public:
    ReimannSolver();
    double4 operator()(const double4& l, const double4& r) const;
};

}

}

#endif // _REIMANN_SOLVER_H_
