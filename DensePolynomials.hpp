#ifndef POLYNOMIALSHPC_DENSEPOLYNOMIALS_HPP
#define POLYNOMIALSHPC_DENSEPOLYNOMIALS_HPP

#include "utils.hpp"
#include <algorithm>
#include <complex>
#include <iostream>
#include <random>
#include <ranges>
#include <type_traits>
#include <vector>

namespace Polynomial {
    template<typename T>
    class Dense;

    template<typename T1, typename T2>
    inline Dense<decltype(T1() * T2())> operator+(const Dense<T1> &,
                                                  const Dense<T2> &);

    template<typename T1, typename T2>
    inline Dense<decltype(T1() * T2())> operator-(const Dense<T1> &,
                                                  const Dense<T2> &);

    template<typename T1, typename T2>
    inline Dense<decltype(T1() * T2())> operator*(const Dense<T1> &,
                                                  const Dense<T2> &);

    template<typename T1, typename T2>
    inline Dense<decltype(T1() * T2())> operator/(const Dense<T1> &,
                                                  const Dense<T2> &);

    template<typename T1, typename T2>
    inline Dense<decltype(T1() * T2())> operator/(const Dense<T1> &, T2);

    template<typename T1, typename T2>
    inline Dense<decltype(T1() * T2())> operator%(const Dense<T1> &,
                                                  const Dense<T2> &);

    template<typename T1, typename T2>
    inline Dense<decltype(T1() * T2())> operator+=(Dense<T1> &, const Dense<T2> &);

    template<typename T1, typename T2>
    inline Dense<decltype(T1() * T2())> operator-=(Dense<T1> &, const Dense<T2> &);

    template<typename T1, typename T2>
    inline Dense<decltype(T1() * T2())> operator*=(Dense<T1> &, const Dense<T2> &);

// template<typename T1, typename T2>
// bool operator==(const Polynomial<T1> &, const Polynomial<T2> &);
//
// template<typename T>
// Polynomial<T> operator/=(const Polynomial<T> &, const Polynomial<T> &);
//
// template<typename T>
// Polynomial<T> operator/=(const Polynomial<T> &, const T);
//
// template<typename T>
// Polynomial<T> operator%=(const Polynomial<T> &, const Polynomial<T> &);

    template<typename T>
    inline std::ostream &operator<<(std::ostream &, const Dense<T> &);

    template<typename T>
    inline Dense<T> derivative(Dense<T> P, int k = 1);

    template<typename T1, typename T2>
    inline Dense<decltype(T1() * T2())> fftmultiply(const Dense<T1> &a,
                                                    const Dense<T2> &b);

    template<typename T>
    class Dense {
    private:
        std::vector<T> coefficients;
        int n;

        //  extend method facilitates certain operations
        [[maybe_unused]] [[nodiscard]] Dense extend(int) const;

        //  adjust removes trailing zeros
        void adjust();

        template<typename T1, typename T2>
        static void multiplyRec(const std::vector<T1> &A, size_t startA, size_t endA,
                                const std::vector<T2> &B, size_t startB, size_t endB,
                                std::vector<decltype(T1() * T2())> &result,
                                size_t startResult) {
            if (endA - startA <= 64 || endB - startB <= 64) {
                // Base case: standard polynomial multiplication
                for (size_t i = startA; i < endA; ++i) {
                    for (size_t j = startB; j < endB; ++j) {
                        result[startResult + i + j - startA - startB] += A[i] * B[j];
                    }
                }
                return;
            }

            size_t midA = std::midpoint(startA, endA);
            size_t midB = std::midpoint(startB, endB);

            // Recursive computations
            multiplyRec(A, startA, midA, B, startB, midB, result, startResult);
            multiplyRec(A, midA, endA, B, midB, endB, result,
                        startResult + midA + midB);
            addAndMultiply(A, startA, midA, midA, endA, B, startB, midB, midB, endB,
                           result, startResult + midA);
        }

        // Helper function to compute (A0 + A1) * (B0 + B1)
        template<typename T1, typename T2>
        static void addAndMultiply(const std::vector<T1> &A, size_t startA0,
                                   size_t endA0, size_t startA1, size_t endA1,
                                   const std::vector<T2> &B, size_t startB0,
                                   size_t endB0, size_t startB1, size_t endB1,
                                   std::vector<decltype(T1() * T2())> &result,
                                   size_t startResult) {
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

            size_t m = n / 3; // Divide size by three

            // Split A and B into three parts each
            size_t midA1 = startA + m;
            size_t midA2 = startA + 2 * m;
            size_t midB1 = startB + m;
            size_t midB2 = startB + 2 * m;

            // Compute intermediate products using Cantor's algorithm
            cantorRec(A, startA, midA1, B, startB, midB1, result, startResult);
            cantorRec(A, midA2, endA, B, midB2, endB, result, startResult + 2 * m);

            std::vector<T> tempA(m);
            std::vector<T> tempB(m);

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

            // Add/subtract intermediate results into the final result, adjusting
            // coefficients
            for (size_t i = 0; i < 2 * m; ++i) {
                result[startResult + m + i] +=
                        result[startResult + i] - result[startResult + 2 * m + i];
            }
        }

    public:
        //    Constructors
        //    Default polynomial is 0, degree=-1 for convenience
        Dense() : coefficients(0), n(-1) {}

        //    Single coefficient polynomial
        explicit Dense(const T &a, const int &m = 0) : coefficients(m + 1, 0), n(m) {
            coefficients[n] = a;
            adjust();
        }

        explicit Dense(const T &&a, const int &m = 0) : coefficients(m + 1, 0), n(m) {
            coefficients[n] = std::move(a);
            adjust();
        }

        explicit Dense(const std::vector<T> &coeffs)
                : coefficients(coeffs), n(coeffs.size() - 1) {
            adjust();
        }

        explicit Dense(const std::vector<T> &&coeffs)
                : coefficients(std::move(coeffs)), n(coeffs.size() - 1) {
            adjust();
        }

        explicit Dense(const std::vector<T> &roots, bool fromRoots);

        //  Copy/move assignments
        Dense<T> &operator=(const Dense<T> &) = default;

        //    Polynomial<T> &operator=(Polynomial<T> &&) = default;

        //  Destructor
        ~Dense() = default;

        //    iterators
        [[nodiscard]] auto begin() const { return coefficients.begin(); }

        [[nodiscard]] auto end() const { return coefficients.end(); }

        //  Basic accessors/mutators
        [[nodiscard]] int degree() const { return n; }

        [[nodiscard]] T dominant() const { return coefficients[n]; }

        //  return image of U by P
        template<typename U>
        T operator()(const U &) const;

        template<typename U>
        std::vector<T> operator()(const std::vector<U> &) const;

        [[nodiscard]] std::vector<T>
        RawMultipointEval(const std::vector<T> &points) const;

        [[nodiscard]] std::vector<T>
        TreeMultipointEval(const std::vector<T> &points) const;

        void derivative(int k = 1);

        [[nodiscard]] T &dominant() { return coefficients[n]; }

        [[nodiscard]] bool is_sparse() const;

        T operator[](int i) const { return (i >= 0 && i <= n) ? coefficients[i] : 0; }

        //  Basic operators as friend methods
        template<typename T1, typename T2>
        friend Dense<decltype(T1() * T2())> operator+(const Dense<T1> &,
                                                      const Dense<T2> &);

        template<typename T1, typename T2>
        friend Dense<decltype(T1() * T2())> operator-(const Dense<T1> &,
                                                      const Dense<T2> &);

        template<typename T1, typename T2>
        friend Dense<decltype(T1() * T2())> operator*(const Dense<T1> &,
                                                      const Dense<T2> &);

        template<typename T1, typename T2>
        friend Dense<decltype(T1() * T2())> operator/(const Dense<T1> &,
                                                      const Dense<T2> &);

        template<typename T1, typename T2>
        friend Dense<decltype(T1() * T2())> operator/(const Dense<T1> &, T2);

        template<typename T1, typename T2>
        friend Dense<decltype(T1() * T2())> operator%(const Dense<T1> &,
                                                      const Dense<T2> &);

        template<typename T1, typename T2>
        friend Dense<decltype(T1() * T2())> operator+=(Dense<T1> &,
                                                       const Dense<T2> &);

        template<typename T1, typename T2>
        friend Dense<decltype(T1() * T2())> operator-=(Dense<T1> &,
                                                       const Dense<T2> &);

        template<typename T1, typename T2>
        friend Dense<decltype(T1() * T2())> operator*=(Dense<T1> &,
                                                       const Dense<T2> &);

        template<typename T1, typename T2>
        friend Dense<decltype(T1() * T2())> operator/=(Dense<T1> &,
                                                       const Dense<T2> &);

        template<typename T1, typename T2>
        friend Dense<decltype(T1() * T2())> operator%=(Dense<T1> &,
                                                       const Dense<T2> &);

        template<typename T1, typename T2>
        friend bool operator==(const Dense<T1> &, const Dense<T2> &);

        friend std::ostream &operator
        <<<>(std::ostream &, const Dense<T> &);

        template<typename T1, typename T2>
        friend Dense<decltype(T1() * T2())> fftmultiply(const Dense<T1> &a,
                                                        const Dense<T2> &b);

        //    template<typename T1, typename T2>
        //    friend Polynomial<decltype(T1() * T2())> karatsuba(const Polynomial<T1>
        //    &A, const Polynomial<T2> &B) {
        //        std::vector<decltype(T1() * T2())> result(A.coefficients.size() +
        //        B.coefficients.size() - 1, 0); multiplyRec(A.coefficients, 0,
        //        A.coefficients.size(),
        //                    B.coefficients, 0, B.coefficients.size(),
        //                    result, 0);
        //        return Polynomial<decltype(T1() * T2())>(result);
        //    }

        //    template<typename T1, typename T2>
        //    friend Polynomial<decltype(T1() * T2())> cantor(const Polynomial<T1> &a,
        //    const Polynomial<T2> &b){
        //                std::vector<decltype(T1() * T2())>
        //                result(a.coefficients.size() + b.coefficients.size() - 1,
        //                0);
        //        cantorRec(a.coefficients, 0, a.coefficients.size(), b.coefficients,
        //        0, b.coefficients.size(), result, 0); return
        //        Polynomial<decltype(T1() * T2())>(result);
        //    };
    };

    template<typename coeffType, typename precision = double>
    std::vector<precision> solveRealRoots(const Dense<coeffType> &P) {
        /*
         * TODO: check GSL implementation
         */
        std::vector<precision> roots;

        if (P.degree() == 1) {
            roots.push_back(static_cast<precision>(-P[0]) /
                            static_cast<precision>(P[1]));
        }

        if (P.degree() == 2) {
            auto delta = P[1] * P[1] - 4 * P[2] * P[0];
            if (delta > 0) {
                roots.push_back((-static_cast<precision>(P[1]) -
                                 sqrt(static_cast<precision>(delta))) /
                                (2 * static_cast<precision>(P[2])));
                roots.push_back((-static_cast<precision>(P[1]) +
                                 sqrt(static_cast<precision>(delta))) /
                                (2 * static_cast<precision>(P[2])));
            }
        }

        if (P.degree() == 3) {
            // delta=18abcd - 4b^3d + b^2c^2 -4ac^3 -27a^2d^2
            // a=P[3], b=P[2], c=P[1], d=P[0]
            auto delta = 18 * P[0] * P[1] * P[2] * P[3] -
                         4 * P[2] * P[2] * P[2] * P[0] + P[2] * P[2] * P[1] * P[1] -
                         4 * P[3] * P[1] * P[1] * P[1] - 27 * P[3] * P[3] * P[0] * P[0];
            if (delta == 0) {
                roots.push_back(static_cast<precision>(
                                        (9 * P[3] * P[0] - P[1] * P[2]) /
                                        static_cast<precision>(2 * P[2] * P[2] - 6 * P[3] * P[1])));
                roots.push_back(static_cast<precision>(
                                        (9 * P[3] * P[0] - P[1] * P[2]) /
                                        static_cast<precision>(2 * P[2] * P[2] - 6 * P[3] * P[1])));
                roots.push_back(static_cast<precision>(
                                        (4 * P[3] * P[2] * P[1] - 9 * P[3] * P[3] * P[0] -
                                         P[2] * P[2] * P[2]) /
                                        static_cast<precision>(P[3] * P[2] * P[2] - 3 * P[3] * P[3] * P[1])));
            }
            /*TODO*/
            if (delta > 0) {
            }
        }

        return roots;
    }

    template<typename coefficientType, typename precision = double>
    std::vector<std::complex<precision>>
    solveComplexRoots(const Dense<coefficientType> &polynomial) {
        // Implementation for solving complex roots with user-specified precision
        std::vector<std::complex<precision>> roots;
        return roots;
    }

    template<typename coefficientType, typename precision = double,
            typename returnType = precision>
    std::vector<returnType> solveRoots(const Dense<coefficientType> &polynomial) {
        // Deduce the return type based on whether the polynomial has integer or
        // floating-point roots
        if constexpr (std::is_floating_point_v<returnType> ||
                      std::is_integral_v<returnType>) {
            return solveRealRoots<coefficientType, returnType>(
                    polynomial); // Solve for real roots with integer coefficients
        } else {
            return solveComplexRoots<coefficientType, returnType>(
                    polynomial); // Solve for complex roots with floating-point coefficients
        }
    }

    template<typename T>
    Dense<T> interpolate(const std::vector<std::pair<T, T>> &points) {
        /*
         * Points are expected to be in increasing order (for the Xs) and distinct.
         * This implementation does not check those conditions
         */
        if (points.empty()) {
            return Dense<T>();
        }
        if (points.size() == 1) {
            return Dense<T>(points[0].second);
        }
        auto m = points.size();
        Dense<T> ans;
        for (long unsigned int i = 0; i < m; ++i) {
            Dense<T> temp({1});
            for (long unsigned int j = 0; j < m; ++j) {
                if (j != i) {
                    temp *= Dense<T>(std::vector<T>{
                            -points[j].first / (points[i].first - points[j].first),
                            1 / (points[i].first - points[j].first)});
                }
            }
            ans += Dense<T>({points[i].second}) * temp;
        }
        return ans;
    }

    template<typename T>
    Dense<T>::Dense(const std::vector<T> &roots, const bool fromRoots) {
        /*
         * creating a polynomial from roots if bool parameter set to True
         */
        if (fromRoots) {
            *this = Dense<T>({1});
            for (const auto &root: roots) {
                *this = *this * Dense<T>(std::vector<T>{-root, 1});
            }
        } else {
            *this = Dense<T>(roots);
        }
        n = coefficients.size() - 1;
    }

    template<typename T>
    bool Dense<T>::is_sparse() const {
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
    T Dense<T>::operator()(const U &x) const {

        /*
         * Using Horner's method for single point evaluation O(n)
         * Possibility to add different methods will be added for sparse cases
         */

        if (n > 0) {
            auto horner = [&x](T s, T a) { return x * s + a; };
            T res = std::accumulate(std::next(coefficients.rbegin()),
                                    coefficients.rend(), coefficients[n], horner);
            return res;
        } else if (n == 0)
            return coefficients[0];
        else
            return 0;
    }

    template<typename T>
    template<typename U>
    std::vector<T> Dense<T>::operator()(const std::vector<U> &v) const {
        // O(n^2)
        //  will try sub quadratic algo with FFT, also need to switch to templated
        //  input vect type
        return RawMultipointEval(v);
    }

    template<typename T>
    std::vector<T> Dense<T>::RawMultipointEval(const std::vector<T> &points) const {

        /*
         * O(n^2) implementation, Horner's method repeated for each point
         */

        int m = points.size();
        std::vector<T> ans(m);
        std::ranges::transform(points, ans.begin(),
                               [this](const T &x) { return this->operator()(x); });
        return ans;
    }

    template<typename T>
    void buildTree(std::vector<Dense<T>> &tree, const std::vector<T> &points,
                   size_t k, size_t idx = 0) {

        /*
         * Helper function generating tree structure for fast multipoint evaluation
         */

        if (idx >= tree.size() / 2) {

            if (size_t pointIdx = idx - tree.size() / 2 < points.size()) {
                tree[idx] = Dense<T>{-points[pointIdx], 1}; // M_0,j = X - u_j
            }
            return;
        }

        buildTree(tree, points, k, 2 * idx + 1);
        buildTree(tree, points, k, 2 * idx + 2);
        if (idx > 0) {
            tree[idx] =
                    tree[2 * idx + 1] * tree[2 * idx + 2]; // M_i,j = M_i-1,2j * M_i-1,2j+1
        }
    }

    template<typename T>
    std::vector<Dense<T>> evaluate(const Dense<T> &f,
                                   const std::vector<Dense<T>> &tree, size_t k,
                                   size_t idx = 0) {

        /*
         * Helper function recursively going through the tree for the fast multipoint
         * evaluation
         */

        std::vector<Dense<T>> result;

        if (idx >= tree.size() / 2) {
            if (idx < tree.size()) {
                result.push_back(f); // Base case
            }
            return result;
        }

        Dense<T> r0 = f % tree[2 * idx + 1]; // f rem M_k-1,0
        Dense<T> r1 = f % tree[2 * idx + 2]; // f rem M_k-1,1

        auto res0 = evaluate(r0, tree, k, 2 * idx + 1);
        auto res1 = evaluate(r1, tree, k, 2 * idx + 2);

        result.insert(result.end(), res0.begin(), res0.end());
        result.insert(result.end(), res1.begin(), res1.end());

        return result;
    }

    template<typename T>
    std::vector<T>
    Dense<T>::TreeMultipointEval(const std::vector<T> &points) const {

        /*
         * Tree based with sub products algorithm for O(nlog(n)) multipoint evaluation
         */

        if (points.empty()) {
            return {};
        }
        if (points.size() == 1) {
            return {this->operator()(points[0])};
        }
        auto m = points.size();
        size_t TreeSize = (1 << (m + 1)) - 1;
        std::vector<Dense<T>> tree(TreeSize);
        buildTree(tree, points, m);
        std::vector<T> ans;

        return ans;
    }

// template<typename T>
// std::vector<T> Polynomial<T>::RawMultipointEval(const std::vector<T> &points)
// const {
//     /*
//      * Algorithm inspired by Modern Computer Algebra - Joachim von zur
//      Gathen, Jürgen Gerhard  (2013)
//      * this is much slower than using Horner n times unfortunately...
//      */
//
//     int m = points.size();
//     if (m <= 10) {
//         std::vector<T> ans(m);
//         std::transform(points.begin(), points.end(), ans.begin(),
//         [this](const T &x) { return this->operator()(x); }); return ans;
//     }
//
//     //Splitting in half
//     std::vector<T> pointsLeft(points.begin(), points.begin() + m / 2);
//     std::vector<T> pointsRight(points.begin() + m / 2, points.end());
//
//     // Construct polynomials from roots
//     Polynomial<T> Q1(pointsLeft, true);
//     Polynomial<T> Q2(pointsRight, true);
//
//     // Compute remainders
//     Polynomial<T> remainderQ1 = (*this) % Q1;
//     Polynomial<T> remainderQ2 = (*this) % Q2;
//
//     // Merge results
//     std::vector<T> results = Q1.RawMultipointEval(pointsLeft);
//     std::vector<T> rightResults = Q2.RawMultipointEval(pointsRight);
//     results.insert(results.end(), rightResults.begin(), rightResults.end());
//     return results;
// }

    inline Dense<int> generateRandomIntPolynomial(const int &n, int inf, int sup) {
        std::random_device rd;
        std::mt19937 G(rd());
        std::uniform_int_distribution unif(inf, sup);
        auto gen = [&G, &unif]() { return unif(G); };
        std::vector<int> v(n + 1);
        std::ranges::generate(v, gen);
        return Dense<int>(v);
    }

    template<typename T, typename Distribution, typename Generator>
    inline Dense<T> generateRandomPolynomial(int degree, Distribution &distribution,
                                             Generator &generator) {
        std::random_device rd;
        auto gen = [&generator, &distribution]() { return distribution(generator); };
        std::vector<T> v(degree + 1);
        std::ranges::generate(v, gen);
        return Dense<T>(v);
    }

    template<typename T>
    void Dense<T>::adjust() {
        /*
         * The while loop may not be ideal
         */
        while ((n >= 0) && (is_zero(coefficients[n]))) {
            n--;
        }
        coefficients.resize(n + 1);
    }

    template<typename T>
    [[maybe_unused]] Dense<T> Dense<T>::extend(int m) const {
        if (m > n) {
            coefficients.reserve(m + 1);
        }
        return this;
    }

    template<typename T>
    Dense<T> derivative(Dense<T> P, int k) {
        P.derivative(k);
        return P;
    }

    template<typename T>
    void Dense<T>::derivative(int k) {
        if (k == 0) [[unlikely]] {
            return;
        }
        if (n < k) {
            *this = Dense<T>();
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
    Dense<decltype(T1() * T2())> operator+(const Dense<T1> &p, const Dense<T2> &q) {
        auto polynomial_ref_pair =
                std::minmax(p, q, [](const Dense<T1> &u, const Dense<T2> &v) {
                    return u.degree() < v.degree();
                });

        Dense<decltype(T1() * T2())> result(
                polynomial_ref_pair.second); // single copy of the largest

        std::ranges::transform(polynomial_ref_pair.first.coefficients,
                               result.coefficients, result.coefficients.begin(),
                               std::plus<decltype(T1() * T2())>());

        result.adjust();

        return result;
    }

    template<typename T1, typename T2>
    Dense<decltype(T1() * T2())> operator-(const Dense<T1> &p, const Dense<T2> &q) {
        auto polynomial_ref_pair =
                std::minmax(p, q, [](const Dense<T1> &u, const Dense<T2> &v) {
                    return u.degree() < v.degree();
                });
        Dense<decltype(T1() * T2())> result(polynomial_ref_pair.second);
        if (&polynomial_ref_pair.first == &p) {
            std::ranges::transform(polynomial_ref_pair.first.coefficients,
                                   result.coefficients, result.coefficients.begin(),
                                   std::minus<decltype(T1() * T2())>());
        } else {
            std::ranges::transform(
                    result.coefficients, polynomial_ref_pair.first.coefficients,
                    result.coefficients.begin(), std::minus<decltype(T1() * T2())>());
        }
        result.adjust();
        result.dominant() = polynomial_ref_pair.second.dominant();
        return result;

        /*
         * The below alternative is slightly slower (but much more readable...)
         * as std::transform is used twice (one time in the plus operation)
         */
        //    Polynomial<decltype(T1() * T2())> result(q);
        //    std::transform(result.coefficients.begin(), result.coefficients.end(),
        //    result.coefficients.begin(),
        //                   std::negate<decltype(T1() * T2())>());
        //    return p + result;
    }

    template<typename T1, typename T2>
    Dense<decltype(T1() * T2())> operator*(const Dense<T1> &p, const Dense<T2> &q) {
        if (p.n < 0 || q.n < 0) [[unlikely]] {
            return Dense<decltype(T1() * T2())>();
        }
        int m = p.n + q.n;
        /*
         * FFT multiplication is accurate for ints only, and become (much) faster for
         * output degree>200
         */
        if (std::is_integral_v<decltype(T1() * T2())> && m > 200) {
            return fftmultiply(p, q);
        } else {
            std::vector<decltype(T1() * T2())> coeffs(m + 1);
            for (int i = 0; i <= p.degree(); ++i) {
                for (int j = 0; j <= q.degree(); ++j) {
                    coeffs[i + j] += p[i] * q[j];
                }
            }
            return Dense<decltype(T1() * T2())>(std::move(coeffs));
        }
    }

    template<typename T>
    std::pair<Dense<T>, Dense<T>> euclid_div(const Dense<T> &a, const Dense<T> &b) {
        Dense<T> q;
        Dense<T> r(a);

        while (r.degree() >= b.degree()) {
            int n = r.degree() - b.degree();
            Dense<T> q1(r.dominant() / b.dominant(), n);
            r = r - b * q1;
            q = q + q1;
        }
        std::pair<Dense<T>, Dense<T>> ans;
        ans.first = std::move(q);
        ans.second = std::move(r);
        return ans;
    }

    template<typename T1, typename T2>
    Dense<decltype(T1() * T2())> operator/(const Dense<T1> &p, const Dense<T2> &q) {
        return euclid_div(p, q).first;
    }

    template<typename T1, typename T2>
    Dense<decltype(T1() * T2())> operator/(const Dense<T1> &p, const T2 q) {
        std::vector<decltype(T1() * T2())> ans(p.coefficients.size());
        std::ranges::transform(p.coefficients, ans.begin(),
                               [q](const T1 x) { return x / q; });
        return Dense(ans);
    }

    template<typename T1, typename T2>
    Dense<decltype(T1() * T2())> operator%(const Dense<T1> &p, const Dense<T2> &q) {
        return euclid_div(p, q).second;
    }

    template<typename T1, typename T2>
    Dense<decltype(T1() * T2())> operator+=(Dense<T1> &lhs, const Dense<T2> &rhs) {
        lhs = lhs + rhs;
        return lhs;
    }

    template<typename T1, typename T2>
    Dense<decltype(T1() * T2())> operator-=(Dense<T1> &lhs, const Dense<T2> &rhs) {
        lhs = lhs - rhs;
        return lhs;
    }

    template<typename T1, typename T2>
    Dense<decltype(T1() * T2())> operator*=(Dense<T1> &lhs, const Dense<T2> &rhs) {
        lhs = lhs * rhs;
        return lhs;
    }

    template<typename T1, typename T2>
    Dense<decltype(T1() * T2())> operator/=(Dense<T1> &lhs, const Dense<T2> &rhs) {
        lhs = lhs / rhs;
        return lhs;
    }

    template<typename T1, typename T2>
    Dense<decltype(T1() * T2())> operator%=(Dense<T1> &lhs, const Dense<T2> &rhs) {
        lhs = lhs % rhs;
        return lhs;
    }

    template<typename T1, typename T2>
    bool operator==(const Dense<T1> &P, const Dense<T2> &Q) {
        if (P.degree() != Q.degree()) {
            return false;
        }
        for (int i = 0; i < P.degree() + 1; i++) {
            if (P[i] != Q[i]) {
                return false;
            }
        }
        return true;
    }

    template<typename T>
    std::ostream &operator<<(std::ostream &out, const Dense<T> &p) {
        char var;
        if constexpr (is_complex<T>) {
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

                    if (i != p.degree()) {
                        if (should_add_plus(p.coefficients[i])) {
                            out << " +";
                        }
                        if (should_add_plus(-p.coefficients[i])) {
                            out << " - ";
                        } else {
                            out << " ";
                        }
                    }

                    if (!is_one(p.coefficients[i]) && !is_one(-p.coefficients[i])) {
                        out << std::abs(p.coefficients[i]);
                    }
                    if (i > 1) {
                        out << var << "^" << i;
                    }
                    if (i == 1) {
                        out << var;
                    }
                }
            }
            if (!is_zero(p.coefficients[0])) {
                if (p.degree() != 0) {
                    if (should_add_plus(p.coefficients[0])) {
                        out << " +";
                    }
                    if (should_add_plus(-p.coefficients[0])) {
                        out << " - ";
                    } else {
                        out << " ";
                    }
                }
                out << std::abs(p.coefficients[0]);
            }
        }
        return out;
    }

    template<typename T1, typename T2>
    Dense<decltype(T1() * T2())> fftmultiply(const Dense<T1> &a,
                                             const Dense<T2> &b) {
        if (a.n < 0 || b.n < 0) {
            return Dense<decltype(T1() * T2())>();
        } else {
            std::vector<std::complex<double>> fa(a.coefficients.begin(),
                                                 a.coefficients.end());
            std::vector<std::complex<double>> fb(b.coefficients.begin(),
                                                 b.coefficients.end());
            long unsigned int n = 1;
            while (n < a.coefficients.size() + b.coefficients.size()) {
                n <<= 1;
            }
            fa.resize(n);
            fb.resize(n);

            fft(fa, false);
            fft(fb, false);
            for (long unsigned int i = 0; i < n; i++) {
                fa[i] *= fb[i];
            }
            fft(fa, true);

            Dense<decltype(T1() * T2())> result;
            result.coefficients.reserve(n);
            for (long unsigned int i = 0; i < n; i++) {
                result.coefficients[i] = round(fa[i].real());
            }
            result.n = a.n + b.n;
            return result;
        }
    }

} // namespace Polynomial

#endif // POLYNOMIALSHPC_DENSEPOLYNOMIALS_HPP