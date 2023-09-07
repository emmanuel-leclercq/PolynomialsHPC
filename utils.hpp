
#pragma once

#include <complex>

template<typename T>
bool is_zero(T a) {
    return (static_cast<int>(a) == 0);
}

template<typename T>
bool is_one(T a) {
    return (a == 1);
}

//specialization for doubles: very small coefficients is basically zero
template<>
bool is_zero(double a) {
    return (std::abs(a) < 1.e-5);
}

template<>
bool is_one<double>(double a) {
    return (std::abs(a - 1.) < 1.e-5);
}

// overload for complex numbers
template<typename T>
bool is_one(std::complex<T> c) {
    return is_one<T>(c.real()) && is_zero<T>(c.imag());
}

template<typename T>
bool is_zero(std::complex<T> c) {
    return is_zero<T>(c.real()) && is_zero<T>(c.imag());
}

template<typename T>
T factorial(int n) {
    T result = 1;
    for (int i = 1; i <= n; i++) {
        result *= i;
    }
    return result;
}

template<typename T>
T rangeProduct(int n, int k) {
    auto x = std::minmax(n, k);

    T result = 1;
    for (int i = x.first; i <= x.second; i++) {
        result *= i;
    }
    return result;
}

template<typename T>
bool should_add_plus(const T &x) {
    return x >= 0;
}

template<typename T>
bool should_add_plus(const std::complex<T> &z) {
    return true;
}