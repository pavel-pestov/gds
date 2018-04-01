#ifndef _EULER_CORRECTOR_H_
#define _EULER_CORRECTOR_H_

#include "types.h"

namespace gds {

namespace math {

class EulerCorrector
{
public:
    EulerCorrector();
    double4 operator()(const double4& l, const double4& x, const double4& r, const double tdh);
    double operator()(const double l, const double x, const double r, const double v, const double tdh);
};

}

}

#endif // _EULER_CORRECTOR_H_
