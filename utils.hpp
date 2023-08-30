
#pragma once

#include <complex>

template <typename T>
bool is_zero(T a)
{
    return (static_cast<int>(a) == 0);
}

template <typename T>
bool is_one(T a)
{
    return (a == 1);
}

//specialization for doubles: very small coefficients is basically zero
template <>
bool is_zero(double a)
{
    return (std::abs(a) < 1.e-5);
}

template <>
bool is_one<double>(double a)
{
    return (std::abs(a - 1.) < 1.e-5);
}

// overload for complex numbers
template <typename T>
bool is_one(std::complex<T> c)
{
    return is_one<T>(c.real()) && is_zero<T>(c.imag());
}

template <typename T>
bool is_zero(std::complex<T> c)
{
    return is_zero<T>(c.real()) && is_zero<T>(c.imag());
}