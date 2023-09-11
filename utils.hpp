
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

class Timer {
    std::chrono::time_point<std::chrono::steady_clock> timePoint;
    size_t value;
public:
    void start() { timePoint = std::chrono::steady_clock::now(); }

    void finish() {
        auto curr = std::chrono::steady_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(curr - timePoint);
        value = elapsed.count();
    }

    size_t operator()() const { return value; }
};

int reverse(int num, int lg_n) {
    int res = 0;
    for (int i = 0; i < lg_n; i++) {
        if (num & (1 << i))
            res |= 1 << (lg_n - 1 - i);
    }
    return res;
}

void fft(std::vector<std::complex<double>> &a, bool invert) {
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