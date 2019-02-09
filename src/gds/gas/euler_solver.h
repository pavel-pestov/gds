#ifndef _EULER_SOLVER_H_
#define _EULER_SOLVER_H_

#include "types.h"

namespace gds {

template<typename T>
vector4<T> euler_solver(const vector4<T>& l, const vector4<T>& x, const vector4<T>& r, const T tdh);

}

#endif // _EULER_SOLVER_H_
