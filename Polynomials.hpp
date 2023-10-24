#ifndef POLYNOMIALSHPC_POLYNOMIALS_HPP
#define POLYNOMIALSHPC_POLYNOMIALS_HPP

#include <vector>
#include <random>
#include <complex>
#include <algorithm>
#include <iostream>
#include "utils.hpp"

template<typename T>
class Polynomial;

template<typename T, typename T1, typename T2>
Polynomial<decltype(T1() * T2())> operator+(const Polynomial<T1> &, const Polynomial<T2> &);

template<typename T, typename T1, typename T2>
Polynomial<decltype(T1() * T2())> operator-(const Polynomial<T1> &, const Polynomial<T2> &);

template<typename T1, typename T2>
Polynomial<decltype(T1() * T2())> operator*(const Polynomial<T1> &, const Polynomial<T2> &);

template<typename T, typename T1, typename T2>
Polynomial<decltype(T1() * T2())> operator/(const Polynomial<T1> &, const Polynomial<T> &);

template<typename T, typename T1, typename T2>
Polynomial<decltype(T1() * T2())> operator%(const Polynomial<T1> &, const Polynomial<T2> &);

template<typename T, typename T1, typename T2>
Polynomial<decltype(T1() * T2())> operator+=(const Polynomial<T1> &, const Polynomial<T2> &);

template<typename T, typename T1, typename T2>
Polynomial<decltype(T1() * T2())> operator-=(const Polynomial<T> &, const Polynomial<T> &);

template<typename T, typename T1, typename T2>
Polynomial<decltype(T1() * T2())> operator*=(const Polynomial<T1> &, const Polynomial<T2> &);

template<typename T1, typename T2>
bool operator==(const Polynomial<T1> &, const Polynomial<T2> &);

template<typename T>
Polynomial<T> operator/=(const Polynomial<T> &, const Polynomial<T> &);

template<typename T>
Polynomial<T> operator%=(const Polynomial<T> &, const Polynomial<T> &);

template<typename T>
std::ostream &operator<<(std::ostream &, const Polynomial<T> &);

template<typename T>
Polynomial<T> derivative(const Polynomial<T> &P, int k = 1);

template<typename T1, typename T2>
Polynomial<decltype(T1() * T2())> fftmultiply(const Polynomial<T1> &a, const Polynomial<T2> &b);

template<typename T>
class Polynomial {
private:
    std::vector<T> coefficients;
    int n;

//  extend method facilitates certain operations
    [[nodiscard]] Polynomial extend(int) const;

//  adjust removes trailing zeros
    void adjust();

    template<typename T1, typename T2>
    static void multiplyRec(const std::vector<T1> &A, size_t startA, size_t endA,
                            const std::vector<T2> &B, size_t startB, size_t endB,
                            std::vector<decltype(T1() * T2())> &result, size_t startResult) {
        if (endA - startA <= 64 || endB - startB <= 64) {
            // Base case: standard polynomial multiplication
            for (size_t i = startA; i < endA; ++i) {
                for (size_t j = startB; j < endB; ++j) {
                    result[startResult + i + j - startA - startB] += A[i] * B[j];
                }
            }
            return;
        }

        size_t midA = (startA + endA) / 2;
        size_t midB = (startB + endB) / 2;

        // Recursive computations
        multiplyRec(A, startA, midA, B, startB, midB, result, startResult);
        multiplyRec(A, midA, endA, B, midB, endB, result, startResult + midA + midB);
        addAndMultiply(A, startA, midA, midA, endA,
                       B, startB, midB, midB, endB,
                       result, startResult + midA);
    }

    // Helper function to compute (A0 + A1) * (B0 + B1)
    template<typename T1, typename T2>
    static void addAndMultiply(const std::vector<T1> &A, size_t startA0, size_t endA0, size_t startA1, size_t endA1,
                               const std::vector<T2> &B, size_t startB0, size_t endB0, size_t startB1, size_t endB1,
                               std::vector<decltype(T1() * T2())> &result, size_t startResult) {
        for (size_t i = startA0; i < endA1; ++i) {
            for (size_t j = startB0; j < endB1; ++j) {
                T aVal = (i < endA0 ? A[i] : A[i - (endA0 - startA0)]);
                T bVal = (j < endB0 ? B[j] : B[j - (endB0 - startB0)]);
                result[startResult + i + j - startA0 - startB0] += aVal * bVal;
            }
        }

        // Subtract the terms that were added twice
        for (size_t i = startA0; i < endA0; ++i) {
            for (size_t j = startB0; j < endB0; ++j) {
                result[startResult + i + j - startA0 - startB0] -= A[i] * B[j];
            }
        }

        for (size_t i = startA1; i < endA1; ++i) {
            for (size_t j = startB1; j < endB1; ++j) {
                result[startResult + i + j - startA0 - startB0] -= A[i] * B[j];
            }
        }
    }

    static void cantorRec(const std::vector<T> &A, size_t startA, size_t endA,
                          const std::vector<T> &B, size_t startB, size_t endB,
                          std::vector<T> &result, size_t startResult) {
        size_t n = endA - startA;

        if (n <= 32) {
            // Base case: standard polynomial multiplication
            for (size_t i = startA; i < endA; ++i) {
                for (size_t j = startB; j < endB; ++j) {
                    result[startResult + i + j - startA - startB] += A[i] * B[j];
                }
            }
            return;
        }

        size_t m = n / 3;  // Divide size by three

        // Split A and B into three parts each
        size_t midA1 = startA + m, midA2 = startA + 2 * m;
        size_t midB1 = startB + m, midB2 = startB + 2 * m;

        // Compute intermediate products using Cantor's algorithm
        cantorRec(A, startA, midA1, B, startB, midB1, result, startResult);
        cantorRec(A, midA2, endA, B, midB2, endB, result, startResult + 2 * m);

        std::vector<T> tempA(m), tempB(m);
        for (size_t i = 0; i < m; ++i) {
            tempA[i] = A[startA + i] + A[midA1 + i] - A[midA2 + i];
            tempB[i] = B[startB + i] + B[midB1 + i] - B[midB2 + i];
        }
        cantorRec(tempA, 0, m, tempB, 0, m, result, startResult + m);

        for (size_t i = 0; i < m; ++i) {
            tempA[i] = A[startA + i] + A[midA1 + i];
            tempB[i] = B[startB + i] + B[midB1 + i];
        }
        cantorRec(tempA, 0, m, tempB, 0, m, result, startResult);

        for (size_t i = 0; i < m; ++i) {
            tempA[i] = A[midA1 + i] - A[midA2 + i];
            tempB[i] = B[midB1 + i] - B[midB2 + i];
        }
        cantorRec(tempA, 0, m, tempB, 0, m, result, startResult + 2 * m);

        // Add/subtract intermediate results into the final result, adjusting coefficients
        for (size_t i = 0; i < 2 * m; ++i) {
            result[startResult + m + i] += result[startResult + i] - result[startResult + 2 * m + i];
        }
    }

public:


//    Constructors
//    Default polynomial is 0, degree=-1 for convenience
    Polynomial() : n(-1), coefficients(0) {}

//    Single coefficient polynomial
    explicit Polynomial(const T &a, const int &m = 0) : n(m), coefficients(m + 1, 0) {
        coefficients[n] = a;
        adjust();
    }

    explicit Polynomial(const T &&a, const int &m = 0) : n(m), coefficients(m + 1, 0) {
        coefficients[n] = std::move(a);
        adjust();
    }

    explicit Polynomial(const std::vector<T> &coeffs) : coefficients(coeffs), n(coeffs.size() - 1) { adjust(); }

    explicit Polynomial(const std::vector<T> &&coeffs) : coefficients(std::move(coeffs)),
                                                         n(coeffs.size() - 1) { adjust(); }

    explicit Polynomial(const std::vector<T> &roots, bool fromRoots);

//  Copy/move assignments
    Polynomial<T> &operator=(const Polynomial<T> &) = default;

//    Polynomial<T> &operator=(Polynomial<T> &&) = default;

//  Destructor
    ~Polynomial() = default;

//  Basic accessors/mutators
    [[nodiscard]] int degree() const { return n; }

    [[nodiscard]] T dominant() const { return coefficients[n]; }

//  return image of U by P
    template<typename U>
    T operator()(const U &) const;

    template<typename U>
    std::vector<T> operator()(const std::vector<U> &) const;

    [[nodiscard]] std::vector<T> multipointEval(const std::vector<T> &points) const;


    void derivative(int k = 1);

    [[nodiscard]] T &dominant() { return coefficients[n]; }

    [[nodiscard]] bool is_sparse() const;

    T operator[](int i) const { return (i >= 0 && i <= n) ? coefficients[i] : 0; }

//  Basic operators as friend methods
    template<typename T1, typename T2>
    friend Polynomial<decltype(T1() * T2())> operator+(const Polynomial<T1> &, const Polynomial<T2> &);

    template<typename T1, typename T2>
    friend Polynomial<decltype(T1() * T2())> operator-(const Polynomial<T1> &, const Polynomial<T2> &);

    template<typename T1, typename T2>
    friend Polynomial<decltype(T1() * T2())>
    operator*(const Polynomial<T1> &, const Polynomial<T2> &);

    template<typename T1, typename T2>
    friend Polynomial<decltype(T1() * T2())> operator/(const Polynomial<T1> &, const Polynomial<T2> &);

    template<typename T1, typename T2>
    friend Polynomial<decltype(T1() * T2())> operator%(const Polynomial<T1> &, const Polynomial<T2> &);

    template<typename T1, typename T2>
    friend Polynomial<decltype(T1() * T2())> operator+=(const Polynomial<T1> &, const Polynomial<T2> &);

    template<typename T1, typename T2>
    friend Polynomial<decltype(T1() * T2())> operator-=(const Polynomial<T1> &, const Polynomial<T2> &);

    template<typename T1, typename T2>
    friend Polynomial<decltype(T1() * T2())> operator*=(const Polynomial<T1> &, const Polynomial<T2> &);

    template<typename T1, typename T2>
    friend Polynomial<decltype(T1() * T2())> operator/=(const Polynomial<T1> &, const Polynomial<T2> &);

    template<typename T1, typename T2>
    friend Polynomial<decltype(T1() * T2())> operator%=(const Polynomial<T1> &, const Polynomial<T2> &);

    template<typename T1, typename T2>
    friend bool operator==(const Polynomial<T1> &, const Polynomial<T2> &);

    friend std::ostream &operator
    <<<>(std::ostream &, const Polynomial<T> &);

    friend Polynomial<T> derivative<T>(const Polynomial<T> &P, int k);

    template<typename T1, typename T2>
    friend Polynomial<decltype(T1() * T2())> fftmultiply(const Polynomial<T1> &a, const Polynomial<T2> &b);


//    template<typename T1, typename T2>
//    friend Polynomial<decltype(T1() * T2())> karatsuba(const Polynomial<T1> &A, const Polynomial<T2> &B) {
//        std::vector<decltype(T1() * T2())> result(A.coefficients.size() + B.coefficients.size() - 1, 0);
//        multiplyRec(A.coefficients, 0, A.coefficients.size(),
//                    B.coefficients, 0, B.coefficients.size(),
//                    result, 0);
//        return Polynomial<decltype(T1() * T2())>(result);
//    }


//    template<typename T1, typename T2>
//    friend Polynomial<decltype(T1() * T2())> cantor(const Polynomial<T1> &a, const Polynomial<T2> &b){
//                std::vector<decltype(T1() * T2())> result(a.coefficients.size() + b.coefficients.size() - 1, 0);
//        cantorRec(a.coefficients, 0, a.coefficients.size(), b.coefficients, 0, b.coefficients.size(), result, 0);
//        return Polynomial<decltype(T1() * T2())>(result);
//    };


    template<typename U>
    Polynomial<T> interpolate(const std::vector<U> &);
};

template<typename T>
template<typename U>
Polynomial<T> Polynomial<T>::interpolate(const std::vector<U> &points) {
    if (points.size() == 0) { return Polynomial<T>(); }
    if (points.size() == 1) { return Polynomial<T>(points[0]); }
    else{

    }
}

template<typename T>
Polynomial<T>::Polynomial(const std::vector<T> &roots, bool fromRoots) {
    if (fromRoots) {
        *this = Polynomial<T>({1});
        for (const auto &root: roots) {
            *this = *this * Polynomial<T>({-root, 1});
        }
    } else { *this = Polynomial<T>(); }
    n = coefficients.size() - 1;
}

template<typename T>
bool Polynomial<T>::is_sparse() const {
    int zeros = 0;
    for (auto c: coefficients) {
        if (is_zero(c)) {
            zeros++;
        }
    }
    return (zeros > n / 2);
}

template<typename T>
template<typename U>
T Polynomial<T>::operator()(const U &x) const {
    /*
     * Using Horner's method for single point evaluation
     * Possibility to add different methods will be added for sparse cases
     */

    if (n > 0) {
        auto horner = [&x](T s, T a) { return x * s + a; };
        T res = std::accumulate(std::next(coefficients.rbegin()), coefficients.rend(), coefficients[n], horner);
        return res;
    } else if (n == 0) return coefficients[0];
    else return 0;
}

template<typename T>
template<typename U>
std::vector<T> Polynomial<T>::operator()(const std::vector<U> &) const {
    return std::vector<T>();
}

template<typename T>
std::vector<T> Polynomial<T>::multipointEval(const std::vector<T> &points) const {
    /*
     * Algorithm inspired by Modern Computer Algebra - Joachim von zur Gathen, Jürgen Gerhard  (2013)
     */

    int m = points.size();
    if (m == 1) {
        return std::vector<T>{this->operator()(points[0])};
    }

    //Splitting in half
    std::vector<T> pointsLeft(points.begin(), points.begin() + m / 2);
    std::vector<T> pointsRight(points.begin() + m / 2, points.end());

    // Construct polynomials from roots
    Polynomial<T> Q1(pointsLeft, true);
    Polynomial<T> Q2(pointsRight, true);

    // Compute remainders
    auto pair = euclid_div(*this, Q1 * Q2);
    Polynomial<T> remainderQ1 = (*this) % Q1;
    Polynomial<T> remainderQ2 = (*this) % Q2;

    // Merge results
    std::vector<T> results = pair.first.multipointEval(pointsLeft);
    std::vector<T> rightResults = pair.second.multipointEval(pointsRight);
    results.insert(results.end(), rightResults.begin(), rightResults.end());
    return results;
}

Polynomial<int> generateRandomIntPolynomial(const int &n, int inf, int sup) {
    std::random_device rd;
    std::mt19937 G(rd());
    std::uniform_int_distribution<int> unif(inf, sup);
    auto gen = [&G, &unif]() { return unif(G); };
    std::vector<int> v(n + 1);
    std::generate(v.begin(), v.end(), gen);
    return Polynomial<int>(v);
}

template<typename T, typename Distribution, typename Generator>
Polynomial<T> generateRandomPolynomial(int n, Distribution &distribution, Generator &generator) {
    std::random_device rd;
    auto gen = [&generator, &distribution]() { return distribution(generator); };
    std::vector<T> v(n + 1);
    std::generate(v.begin(), v.end(), gen);
    return Polynomial<T>(v);
}

template<typename T>
void Polynomial<T>::adjust() {
//    Here we optimize condition order by lowest to highest probability,
//    since the second condition will not be checked if first one is false.
//    The while loop may not be ideal
    while ((is_zero(coefficients[n])) && (n >= 0)) {
        n--;
    }
    coefficients.resize(n + 1);
}

template<typename T>
Polynomial<T> Polynomial<T>::extend(int m) const {
    if (m > n) { coefficients.reserve(m + 1); }
    return this;
}

template<typename T>
Polynomial<T> derivative(const Polynomial<T> &P, int k) {
    if (k == 0) [[unlikely]] { return P; }
    if (P.n < k) { return Polynomial<T>(); }
    auto Q = P;
    if (P.n < k) {
        Q.coefficients.resize(Q.n + 1 - k);
        for (int i = Q.n; i > -1; i--) {
            Q.coefficients[i] /= rangeProduct<T>(i + 1, i - k);
            std::swap(Q.coefficients[i], Q.coefficients[i - k]);
        }
        Q.n -= k;
    }
    if (k > 0 && P.n >= k) {
        Q.coefficients.resize(Q.n + 1 - k);
        T a = factorial<T>(k);
        for (int i = 0; i < Q.n + 1 - k; i++) {
            Q.coefficients[i] = Q.coefficients[i + k] * a;
            a = a * (k + i + 1) / (i + 1);
        }
        Q.n -= k;
    }
    Q.adjust();
    return Q;
}

template<typename T>
void Polynomial<T>::derivative(int k) {
    if (k == 0) [[unlikely]] { return; }
    if (n < k) {
        *this = Polynomial<T>();
        return;
    }
    if (k < 0) {
        coefficients.resize(n + 1 - k);
        for (int i = n; i > -1; i--) {
            coefficients[i] /= rangeProduct<T>(i + 1, i - k);
            std::swap(coefficients[i], coefficients[i - k]);
        }
        n -= k;
    }
    if (k > 0 && n >= k) {
        this->coefficients.resize(n + 1 - k);
        T a = factorial<T>(k);
        for (int i = 0; i < n + 1 - k; i++) {
            coefficients[i] = coefficients[i + k] * a;
            a = a * (k + i + 1) / (i + 1);
        }
        n -= k;
    }
    adjust();
}

template<typename T1, typename T2>
Polynomial<decltype(T1() * T2())> operator+(const Polynomial<T1> &p, const Polynomial<T2> &q) {
    auto polynomial_ref_pair = std::minmax(p, q, [](const Polynomial<T1> &u, const Polynomial<T2> &v) {
        return u.degree() < v.degree();
    });

    Polynomial<decltype(T1() * T2())> result(polynomial_ref_pair.second); // single copy of the largest

    std::transform(
            polynomial_ref_pair.first.coefficients.begin(), polynomial_ref_pair.first.coefficients.end(),
            result.coefficients.begin(),
            result.coefficients.begin(),
            std::plus<decltype(T1() * T2())>());

    result.adjust();

    return result;
}

template<typename T1, typename T2>
Polynomial<decltype(T1() * T2())> operator-(const Polynomial<T1> &p, const Polynomial<T2> &q) {
    auto polynomial_ref_pair = std::minmax(p, q, [](const Polynomial<T1> &u, const Polynomial<T2> &v) {
        return u.degree() < v.degree();
    });
    Polynomial<decltype(T1() * T2())> result(polynomial_ref_pair.second);
    if (&polynomial_ref_pair.first == &p) {
        std::transform(
                polynomial_ref_pair.first.coefficients.begin(), polynomial_ref_pair.first.coefficients.end(),
                result.coefficients.begin(),
                result.coefficients.begin(),
                std::minus<decltype(T1() *
                                    T2())>());
    } else {
        std::transform(result.coefficients.begin(),
                       result.coefficients.end(), polynomial_ref_pair.first.coefficients.begin(),
                       result.coefficients.begin(),
                       std::minus<decltype(T1() *
                                           T2())>());
    }
    result.adjust();
    result.dominant() = polynomial_ref_pair.second.dominant();
    return result;

    /*
     * The below alternative is slightly slower (but much more readable...)
     * as std::transform is used twice (one time in the plus operation)
     */
//    Polynomial<decltype(T1() * T2())> result(q);
//    std::transform(result.coefficients.begin(), result.coefficients.end(), result.coefficients.begin(),
//                   std::negate<decltype(T1() * T2())>());
//    return p + result;
}

template<typename T1, typename T2>
Polynomial<decltype(T1() * T2())> operator*(const Polynomial<T1> &p, const Polynomial<T2> &q) {
    if (p.n < 0 || q.n < 0)[[unlikely]] {
        return Polynomial<decltype(T1() * T2())>();
    }
    int m = p.n + q.n;
    /*
     * FFT multiplication is accurate for ints only, and become (much) faster for output degree>200
     */
    if (std::is_same<decltype(T1() * T2()), int>::value && m > 200) { return fftmultiply(p, q); }
    else {
        std::vector<decltype(T1() * T2())> coeffs(m + 1);
        for (int i = 0; i <= p.degree(); ++i) {
            for (int j = 0; j <= q.degree(); ++j) {
                coeffs[i + j] += p[i] * q[j];
            }
        }
        return Polynomial<decltype(T1() * T2())>(std::move(coeffs));
    }
}

template<typename T>
std::pair<Polynomial<T>, Polynomial<T>> euclid_div(const Polynomial<T> &a, const Polynomial<T> &b) {
    Polynomial<T> q;
    Polynomial<T> r(a);

    while (r.degree() >= b.degree()) {
        int n = r.degree() - b.degree();
        Polynomial<T> q1(r.dominant() / b.dominant(), n);
        r = r - b * q1;
        q = q + q1;
    }
    std::pair<Polynomial<T>, Polynomial<T>> ans;
    ans.first = std::move(q);
    ans.second = std::move(r);
    return ans;
}

template<typename T1, typename T2>
Polynomial<decltype(T1() * T2())> operator/(const Polynomial<T1> &p, const Polynomial<T2> &q) {
    return euclid_div(p, q).first;
}

template<typename T1, typename T2>
Polynomial<decltype(T1() * T2())> operator%(const Polynomial<T1> &p, const Polynomial<T2> &q) {
    return euclid_div(p, q).second;
}

template<typename T1, typename T2>
Polynomial<decltype(T1() * T2())> operator+=(const Polynomial<T1> &lhs, const Polynomial<T2> &rhs) {
    return lhs + rhs;
}

template<typename T1, typename T2>
Polynomial<decltype(T1() * T2())> operator-=(const Polynomial<T1> &lhs, const Polynomial<T2> &rhs) {
    return lhs - rhs;
}

template<typename T1, typename T2>
Polynomial<decltype(T1() * T2())> operator*=(const Polynomial<T1> &lhs, const Polynomial<T2> &rhs) {
    return lhs * rhs;
}

template<typename T1, typename T2>
Polynomial<decltype(T1() * T2())> operator/=(const Polynomial<T1> &lhs, const Polynomial<T2> &rhs) {
    return lhs / rhs;
}

template<typename T1, typename T2>
Polynomial<decltype(T1() * T2())> operator%=(const Polynomial<T1> &lhs, const Polynomial<T2> &rhs) {
    return lhs % rhs;
}

template<typename T1, typename T2>
bool operator==(const Polynomial<T1> &P, const Polynomial<T2> &Q) {
    if (P.degree() != Q.degree()) { return false; }
    for (int i = 0; i < P.degree() + 1; i++) { if (P[i] != Q[i]) { return false; }}
    return true;
}

template<typename T>
std::ostream &operator<<(std::ostream &out, const Polynomial<T> &p) {
    // +/- implementation ?
    char var;
    if constexpr (std::is_same<T, std::complex<float>>::value || std::is_same<T, std::complex<double>>::value ||
                  std::is_same<T, std::complex<long double>>::value) {
        var = 'z';
        for (int i = p.n; i > 1; --i) {
            if (!is_zero(p.coefficients[i])) {
                out << p.coefficients[i] << var << "^" << i << " + ";
            }
        }
        out << p.coefficients[0] << var << "^" << 0;
    } else {
        var = 'x';

        if (p.n < 0) {
            out << "0";
            return out;
        }

        for (int i = p.n; i > 0; --i) {
            if (!is_zero(p.coefficients[i])) {

                if (i != p.degree())[[likely]] {
                    if (should_add_plus(p.coefficients[i])) {
                        out << " +";
                    }
                    if (should_add_plus(-p.coefficients[i])) {
                        out << " - ";
                    } else { out << " "; }
                }

                if (!is_one(p.coefficients[i]) && !is_one(-p.coefficients[i])) {
                    out << std::abs(p.coefficients[i]);
                }
                if (i > 1)[[likely]] { out << var << "^" << i; }
                if (i == 1) { out << var; }
            }
        }
        if (!is_zero(p.coefficients[0])) {
            if (p.degree() != 0) {
                if (should_add_plus(p.coefficients[0])) { out << " +"; }
                if (should_add_plus(-p.coefficients[0])) { out << " - "; }
                else { out << " "; }
            }
            out << std::abs(p.coefficients[0]);
        }
    }
    return out;
}

template<typename T1, typename T2>
Polynomial<decltype(T1() * T2())> fftmultiply(const Polynomial<T1> &a, const Polynomial<T2> &b) {
    if (a.n < 0 || b.n < 0) {
        Polynomial<decltype(T1() * T2())> zero;
        return zero;
    } else {
        std::vector<std::complex<double>> fa(a.coefficients.begin(), a.coefficients.end()), fb(
                b.coefficients.begin(), b.coefficients.end());
        int n = 1;
        while (n < a.coefficients.size() + b.coefficients.size()) { n <<= 1; }
        fa.resize(n);
        fb.resize(n);

        fft(fa, false);
        fft(fb, false);
        for (int i = 0; i < n; i++) { fa[i] *= fb[i]; }
        fft(fa, true);

        Polynomial<decltype(T1() * T2())> result;
        result.coefficients.reserve(n);
        for (int i = 0; i < n; i++) {
            result.coefficients[i] = round(fa[i].real());
        }
        result.n = a.n + b.n;
        return result;
    }
}


#endif //POLYNOMIALSHPC_POLYNOMIALS_HPP