#ifndef _EULER_CORRECTOR_H_
#define _EULER_CORRECTOR_H_

#include "types.h"

namespace gds {

template<typename T>
vector4<T> corrector(const vector4<T>& l, const vector4<T>& x, const vector4<T>& r, const T tdh);

template<typename T>
T corrector(const T l, const T x, const T r, const T v, const T tdh);

}

#endif // _EULER_CORRECTOR_H_
