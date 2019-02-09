#ifndef _REIMANN_SOLVER_H_
#define _REIMANN_SOLVER_H_

#include "types.h"

namespace gds {

template<typename T>
vector4<T> rieman_solver(const vector4<T>& l, const vector4<T>& r);

}

#endif // _REIMANN_SOLVER_H_
