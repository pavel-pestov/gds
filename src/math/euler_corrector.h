#ifndef _EULER_CORRECTOR_H_
#define _EULER_CORRECTOR_H_

#include "types.h"

namespace gds {

namespace math {

class EulerCorrector
{
public:

    enum class Order
    {
        first,
        second_fast,
        second_full,
        third_fast,
        third_full
    };

    EulerCorrector(const Order order);
    double4 operator()(const double4& l, const double4& x, const double4& r, const double tdh);
    double operator()(const double l, const double x, const double r, const double v, const double tdh);
private:
    static double third_order_full(const double dl, const double dr, const double shift);
    static double third_order_fast(const double dl, const double dr, const double shift);
    static double second_order_full(const double dl, const double dr, const double shift);
    static double second_order_fast(const double dl, const double dr, const double shift);
    static double first_order(const double dl, const double dr, const double shift);
    double4 decomposition(const double dro, double dp, double dv, const double ro, const double bc, const double bc2);
    double4 correction(double4 l, double4 x, double4 r, double tdh);
    double4 euler_correction(const double4& l, const double4& x, const double4& r, const double tdh);

    double (*high_order)(double, double, double);
};

}

}

#endif // _EULER_CORRECTOR_H_
