#ifndef _FLUID_CORRECTOR_H_
#define _FLUID_CORRECTOR_H_

#include "types.h"

namespace gds {

namespace math {

class FluidCorrector
{
public:
    FluidCorrector();
    double4 operator()(const double4& l, const double4& x, const double4& r, const float tdh);
};

}

}

#endif // _FLUID_CORRECTOR_H_
