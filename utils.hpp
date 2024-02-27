// utils.hpp
#pragma once

#include <complex>

#include <chrono>


template<typename T>
inline bool is_zero(T a) {
    return (static_cast<int>(a) == 0);
}

template<typename T>
inline bool is_one(T a) {
    return (a == 1);
}

//specialization for doubles: very small coefficients is basically zero
template<>
inline bool is_zero(double a) {
    return (std::abs(a) < 1.e-5);
}

template<>
inline bool is_one<double>(double a) {
    return (std::abs(a - 1.) < 1.e-5);
}

// overload for complex numbers
template<typename T>
inline bool is_one(std::complex<T> c) {
    return is_one<T>(c.real()) && is_zero<T>(c.imag());
}

template<typename T>
inline bool is_zero(std::complex<T> c) {
    return is_zero<T>(c.real()) && is_zero<T>(c.imag());
}

template<typename T>
inline T factorial(int n) {
    T result = 1;
    for (int i = 1; i <= n; i++) {
        result *= i;
    }
    return result;
}

template<typename T>
inline T rangeProduct(int n, int k) {
    auto x = std::minmax(n, k);

    T result = 1;
    for (int i = x.first; i <= x.second; i++) {
        result *= i;
    }
    return result;
}

template<typename T>
inline bool should_add_plus(const T &x) {
    return x >= 0;
}

template<typename T>
inline bool should_add_plus(const std::complex<T> &z) {
    return true;
}

inline int reverse(int num, int lg_n) {
    int res = 0;
    for (int i = 0; i < lg_n; i++) {
        if (num & (1 << i))
            res |= 1 << (lg_n - 1 - i);
    }
    return res;
}

inline void fft(std::vector<std::complex<double>> &a, bool invert) {
    int n = a.size();
    int lg_n = 0;
    while ((1 << lg_n) < n)
        lg_n++;

    for (int i = 0; i < n; i++) {
        if (i < reverse(i, lg_n))
            swap(a[i], a[reverse(i, lg_n)]);
    }

    for (int len = 2; len <= n; len <<= 1) {
        double ang = 2 * M_PI / len * (invert ? -1 : 1);
        std::complex<double> wlen(cos(ang), sin(ang));
        for (int i = 0; i < n; i += len) {
            std::complex<double> w(1);
            for (int j = 0; j < len / 2; j++) {
                std::complex<double> u = a[i + j], v = a[i + j + len / 2] * w;
                a[i + j] = u + v;
                a[i + j + len / 2] = u - v;
                w *= wlen;
            }
        }
    }

    if (invert) {
        for (auto &x: a)
            x /= n;
    }
}

constexpr bool IsPowerOf2(const size_t value) {
    return value && (!(value & (value - 1)));
}

template<typename T>
inline void coutVect(std::ostream &os, const std::vector<T> &v) {
    for (auto x: v) { os << x << " "; }
}

class Timer {
    std::chrono::time_point<std::chrono::steady_clock> timePoint;
    size_t value;
public:
    inline void start() { timePoint = std::chrono::steady_clock::now(); }

    inline void finish() {
        auto curr = std::chrono::steady_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(curr - timePoint);
        value = elapsed.count();
    }

    inline size_t operator()() const { return value; }
};