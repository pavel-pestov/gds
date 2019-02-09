#ifndef _FLUID_CORRECTOR_H_
#define _FLUID_CORRECTOR_H_

#include "types.h"

namespace gds {

template<typename T>
vector4<T> fluid_limiter(const vector4<T>& l, const vector4<T>& x, const vector4<T>& r, const T tdh);

}

#endif // _FLUID_CORRECTOR_H_
