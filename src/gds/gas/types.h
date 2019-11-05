#ifndef _TYPES_H_
#define _TYPES_H_

#include <math.h>
#include <cstring>

namespace gds {

template<typename T> inline T min(const T a, const T b)
{
    return static_cast<T>(a < b) * (a - b) + b;
}

template<typename T> inline T max(const T a, const T b)
{
    return static_cast<T>(a < b) * (b - a) + a;
}

template<typename T> inline T sgn(const T x)
{
    return static_cast<T>(0.0 < x) * 2.0 - 1.0;
}

template<typename T> inline T zero3(const T x)
{
    if (static_cast<T>(1) <= x)
        return static_cast<T>(1);
    if (x <= static_cast<T>(0))
        return static_cast<T>(0);
    return x * (static_cast<T>(2) - x) * static_cast<T>(1.002003004) - static_cast<T>(0.002003004);
}

template<typename T> inline T zero6(const T x)
{
    if (static_cast<T>(1) <= x)
        return static_cast<T>(1);
    if (x <= static_cast<T>(0))
        return static_cast<T>(0);
    return x * (static_cast<T>(2) - x) * static_cast<T>(1.000002003) - static_cast<T>(0.0000002003);
}

template <typename T> inline T sign(const T x) {
    return 1 - ((x < T(0)) << 1);
}

template<typename T>
class Vector4
{
public:
    Vector4()
    {
        _x[0] = _x[1] = _x[2] = _x[3] = 0.0;
    }

    Vector4(const T x)
    {
        _x[0] = _x[1] = _x[2] = _x[3] = x;
    }

    Vector4(const T x0, const T x1, const T x2, const T x3)
    {
        _x[0] = x0;
        _x[1] = x1;
        _x[2] = x2;
        _x[3] = x3;
    }

    Vector4(const T* x)
    {
        memcpy(_x, x, 4 * sizeof(T));
    }

    Vector4& operator=(const T* x)
    {
        memcpy(_x, x, 4 * sizeof(T));
        return *this;
    }

    inline const T* operator*() const
    {
        return _x;
    }

    inline T* operator*()
    {
        return _x;
    }

    inline operator const T*() const
    {
        return _x;
    }

    inline operator T*()
    {
        return _x;
    }

    inline const T operator[](int i) const
    {
        return _x[i];
    }

    inline T& operator[](int i)
    {
        return _x[i];
    }

    inline const T ro() const
    {
        return _x[0];
    }

    inline T& ro()
    {
        return _x[0];
    }

    inline const T p() const
    {
        return _x[1];
    }

    inline T& p()
    {
        return _x[1];
    }

    inline const T v() const
    {
        return _x[2];
    }

    inline T& v()
    {
        return _x[2];
    }

    inline const T g() const
    {
        return _x[3];
    }

    inline T& g()
    {
        return _x[3];
    }

    inline const T i() const
    {
        return ro() * v();
    }

    inline T i(T x)
    {
        _x[2] = x / ro();
        return x;
    }

    inline const T t() const
    {
        return p() / ro();
    }

    inline T t(T x)
    {
        _x[1] = x * ro();
        return x;
    }

    inline const T c() const
    {
        return sqrt(g() * t());
    }

    inline const T ek() const
    {
        return ro() * v() * v() * 0.5;
    }

    inline const T ep() const
    {
        return p() / (g() - 1.0);
    }

    inline T ep(T x)
    {
        _x[1] = x * (g() - 1.0);
        return x;
    }

    inline const T e() const
    {
        return ep() + ek();
    }

    inline T e(T x)
    {
        ep(x - ek());
        return x;
    }

private:
    T _x[4];
};

template<typename T>
class SimpleCell
{
public:
    Vector4<T> gas;
    T h;

    inline const T operator[](int i) const
    {
        return gas[i];
    }

    inline T& operator[](int i)
    {
        return gas[i];
    }
};

}

#endif // _TYPES_H_
