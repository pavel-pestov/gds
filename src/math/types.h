#ifndef _TYPES_H_
#define _TYPES_H_

#include <math.h>

namespace gds {

template<typename T> inline T min(const T a, const T b)
{
    return a < b ? a : b;
}

template<typename T> inline T max(const T a, const T b)
{
    return a < b ? b : a;
}

template<typename T> inline T sgn(const T x)
{
    return x < static_cast<T>(0) ? static_cast<T>(-1) :  static_cast<T>(1);
}

template<typename T> inline T zero3(T x)
{
    if (static_cast<T>(1) <= x)
        return static_cast<T>(1);
    if (x <= static_cast<T>(0))
        return static_cast<T>(0);
    return x * (static_cast<T>(2) - x) * static_cast<T>(1.002003004) - static_cast<T>(0.002003004);
}

template<typename T> inline T zero6(T x)
{
    if (static_cast<T>(1) <= x)
        return static_cast<T>(1);
    if (x <= static_cast<T>(0))
        return static_cast<T>(0);
    return x * (static_cast<T>(2) - x) * static_cast<T>(1.000002003) - static_cast<T>(0.0000002003);
}

template <typename T> T sign(T x) {
    return (T(0) < x) - (x < T(0));
}

class float4
{
public:
    float4()
    {
        _x[0] = _x[1] = _x[2] = _x[3] = 0.0f;
    }

    float4(const float x)
    {
        _x[0] = _x[1] = _x[2] = _x[3] = x;
    }

    float4(const float x0, const float x1, const float x2, const float x3)
    {
        _x[0] = x0;
        _x[1] = x1;
        _x[2] = x2;
        _x[3] = x3;
    }

    float4(const float* x)
    {
        _x[0] = x[0];
        _x[1] = x[1];
        _x[2] = x[2];
        _x[3] = x[3];
    }

    float4(const double* x)
    {
        copy_from(x);
    }

    float4& operator=(const double* x)
    {
        copy_from(x);
        return *this;
    }

    inline const float* operator*() const
    {
        return _x;
    }

    inline float* operator*()
    {
        return _x;
    }

    inline operator const float*() const
    {
        return _x;
    }

    inline operator float*()
    {
        return _x;
    }

    inline const float operator[](int i) const
    {
        return _x[i];
    }

    inline float& operator[](int i)
    {
        return _x[i];
    }

    inline const float ro() const
    {
        return _x[0];
    }

    inline float& ro()
    {
        return _x[0];
    }

    inline const float p() const
    {
        return _x[1];
    }

    inline float& p()
    {
        return _x[1];
    }

    inline const float v() const
    {
        return _x[2];
    }

    inline float& v()
    {
        return _x[2];
    }

    inline const float g() const
    {
        return _x[3];
    }

    inline float& g()
    {
        return _x[3];
    }

    inline const float i() const
    {
        return ro() * v();
    }

    inline float i(float x)
    {
        _x[2] = x / ro();
        return x;
    }

    inline const float t() const
    {
        return p() / ro();
    }

    inline float t(float x)
    {
        _x[1] = x * ro();
        return x;
    }

    inline const float c() const
    {
        return sqrt(g() * t());
    }

    inline const float ek() const
    {
        return ro() * v() * v() * 0.5;
    }

    inline const float ep() const
    {
        return p() / (g() - 1.0);
    }

    inline float ep(float x)
    {
        _x[1] = x * (g() - 1.0);
        return x;
    }

    inline const float e() const
    {
        return ep() + ek();
    }

    inline float e(float x)
    {
        ep(x - ek());
        return x;
    }

private:
    inline void copy_from(const double* x)
    {
        _x[0] = static_cast<float>(x[0]);
        _x[1] = static_cast<float>(x[1]);
        _x[2] = static_cast<float>(x[2]);
        _x[3] = static_cast<float>(x[3]);
    }

    float _x[4];
};

class double4
{
public:
    double4()
    {
        _x[0] = _x[1] = _x[2] = _x[3] = 0.0f;
    }

    double4(const double x)
    {
        _x[0] = _x[1] = _x[2] = _x[3] = x;
    }

    double4(const double x0, const double x1, const double x2, const double x3)
    {
        _x[0] = x0;
        _x[1] = x1;
        _x[2] = x2;
        _x[3] = x3;
    }

    double4(const double* x)
    {
        _x[0] = x[0];
        _x[1] = x[1];
        _x[2] = x[2];
        _x[3] = x[3];
    }

    double4(const float* x)
    {
        copy_from(x);
    }

    double4& operator=(const float* x)
    {
        copy_from(x);
        return *this;
    }

    inline const double* operator*() const
    {
        return _x;
    }

    inline double* operator*()
    {
        return _x;
    }

    inline operator const double*() const
    {
        return _x;
    }

    inline operator double*()
    {
        return _x;
    }

    inline const double operator[](int i) const
    {
        return _x[i];
    }

    inline double& operator[](int i)
    {
        return _x[i];
    }

    inline const double ro() const
    {
        return _x[0];
    }

    inline double& ro()
    {
        return _x[0];
    }

    inline const double p() const
    {
        return _x[1];
    }

    inline double& p()
    {
        return _x[1];
    }

    inline const double v() const
    {
        return _x[2];
    }

    inline double& v()
    {
        return _x[2];
    }

    inline const double g() const
    {
        return _x[3];
    }

    inline double& g()
    {
        return _x[3];
    }

    inline const double i() const
    {
        return ro() * v();
    }

    inline double i(double x)
    {
        _x[2] = x / ro();
        return x;
    }

    inline const double t() const
    {
        return p() / ro();
    }

    inline double t(double x)
    {
        _x[1] = x * ro();
        return x;
    }

    inline const double c() const
    {
        return sqrt(g() * t());
    }

    inline const double ek() const
    {
        return ro() * v() * v() * 0.5;
    }

    inline const double ep() const
    {
        return p() / (g() - 1.0);
    }

    inline double ep(double x)
    {
        _x[1] = x * (g() - 1.0);
        return x;
    }

    inline const double e() const
    {
        return ep() + ek();
    }

    inline double e(double x)
    {
        ep(x - ek());
        return x;
    }

private:
    inline void copy_from(const float* x)
    {
        _x[0] = static_cast<double>(x[0]);
        _x[1] = static_cast<double>(x[1]);
        _x[2] = static_cast<double>(x[2]);
        _x[3] = static_cast<double>(x[3]);
    }

    double _x[4];
};

}

#endif // _TYPES_H_
