#ifndef POLYNOMIALSHPC_POLYNOMIALS_HPP
#define POLYNOMIALSHPC_POLYNOMIALS_HPP

#include <vector>
#include <random>
#include <complex>
#include <algorithm>
#include "utils.hpp"

template<typename T>
class Polynomial;

template<typename T>
Polynomial<T> operator+(const Polynomial<T> &, const Polynomial<T> &);

template<typename T>
Polynomial<T> operator-(const Polynomial<T> &, const Polynomial<T> &);

template<typename T>
Polynomial<T> operator*(const Polynomial<T> &, const Polynomial<T> &);

template<typename T>
Polynomial<T> operator/(const Polynomial<T> &, const Polynomial<T> &);

template<typename T>
Polynomial<T> operator%(const Polynomial<T> &, const Polynomial<T> &);

template<typename T>
Polynomial<T> operator+=(const Polynomial<T> &, const Polynomial<T> &);

template<typename T>
Polynomial<T> operator-=(const Polynomial<T> &, const Polynomial<T> &);

template<typename T>
Polynomial<T> operator*=(const Polynomial<T> &, const Polynomial<T> &);

template<typename T>
Polynomial<T> operator/=(const Polynomial<T> &, const Polynomial<T> &);

template<typename T>
Polynomial<T> operator%=(const Polynomial<T> &, const Polynomial<T> &);

template<typename T>
std::ostream &operator<<(std::ostream &, const Polynomial<T> &);

template<typename T>
class Polynomial {
private:
    std::vector<T> coefficients;
//  Replace degree with template int typename for flexibility?
    int n;

//  extend method facilitates certain operations
    [[nodiscard]] Polynomial extend(int) const;

//  adjust removes trailing zeros
    void adjust();

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

//  Copy/move assignments
    Polynomial<T> &operator=(const Polynomial<T> &) = default;

//    Polynomial<T> &operator=(Polynomial<T> &&) = default;

//  Destructor
    ~Polynomial() = default;


//  return image of U by P
    template<typename U>
    T operator()(const U &);

//  Basic accessors/mutators
    [[nodiscard]] int degree() const { return n; }

    [[nodiscard]] T dominant() const { return coefficients[n]; }

    [[nodiscard]] T &dominant() { return coefficients[n]; }

    [[nodiscard]] bool is_sparse() const;

    T operator[](int i) const { return (i >= 0 && i <= n) ? coefficients[i] : 0; }

//  Basic operators as friend methods
    friend Polynomial<T> operator+<T>(const Polynomial<T> &, const Polynomial<T> &);

    friend Polynomial<T> operator-<T>(const Polynomial<T> &, const Polynomial<T> &);

    friend Polynomial<T> operator*<T>(const Polynomial<T> &, const Polynomial<T> &);

    friend Polynomial<T> operator/<T>(const Polynomial<T> &, const Polynomial<T> &);

    friend Polynomial<T> operator%<T>(const Polynomial<T> &, const Polynomial<T> &);

    friend Polynomial<T> operator+=<T>(const Polynomial<T> &, const Polynomial<T> &);

    friend Polynomial<T> operator-=<T>(const Polynomial<T> &, const Polynomial<T> &);

    friend Polynomial<T> operator*=<T>(const Polynomial<T> &, const Polynomial<T> &);

    friend Polynomial<T> operator/=<T>(const Polynomial<T> &, const Polynomial<T> &);

    friend Polynomial<T> operator%=<T>(const Polynomial<T> &, const Polynomial<T> &);

    friend std::ostream &operator
    <<<>(std::ostream &, const Polynomial<T> &);
};

template<typename T>
bool Polynomial<T>::is_sparse() const {
    int zeros = 0;
    for (auto c: coefficients) {
        if (is_zero(c)) {
            zeros++;
        }
    }
    return (zeros > this->degree() / 2);
}

template<typename T>
template<typename U>
T Polynomial<T>::operator()(const U &x) {
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

/*
 * Will replace uniform int distribution with generic distribution type,
 * and put default (null?) values pour inf and sup.
 * Could also randomize polynomial degree, or leave it to user
 */

Polynomial<double> generateRandomPolynomial(const int &n, double inf, double sup) {
    std::random_device rd;
    std::mt19937 G(rd());
    std::uniform_int_distribution<int> unif(inf, sup);
    auto gen = [&G, &unif]() { return unif(G); };
    std::vector<double> v(n + 1);
    std::generate(v.begin(), v.end(), gen);
    return Polynomial<double>(v);
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
Polynomial<T> operator+(const Polynomial<T> &p, const Polynomial<T> &q) {
    auto polynomial_ref_pair = std::minmax(p, q, [](const Polynomial<T> &u, const Polynomial<T> &v) {
        return u.degree() < v.degree();
    });

    Polynomial<T> result(polynomial_ref_pair.second); // single copy of the largest

    std::transform(
            polynomial_ref_pair.first.coefficients.begin(), polynomial_ref_pair.first.coefficients.end(),
            result.coefficients.begin(),
            result.coefficients.begin(),
            std::plus<T>());

    result.adjust();

    return result;
}

template<typename T>
Polynomial<T> operator-(const Polynomial<T> &p, const Polynomial<T> &q) {
    auto polynomial_ref_pair = std::minmax(p, q, [](const Polynomial<T> &u, const Polynomial<T> &v) {
        return u.degree() < v.degree();
    });
    Polynomial<T> result(polynomial_ref_pair.second); // copie du plus grand dans le résultat.
    if (&polynomial_ref_pair.first == &p) {// pour trouver qui est celui qui porte le - entre le min et le max!
        std::transform(
                polynomial_ref_pair.first.coefficients.begin(), polynomial_ref_pair.first.coefficients.end(),
                result.coefficients.begin(),
                result.coefficients.begin(),
                std::minus<T>()); // l'ordre des arguments permet de ne parcourir que jusqu'au degré du plus petit polynôme.
    } else {
        std::transform(result.coefficients.begin(),
                       result.coefficients.end(), polynomial_ref_pair.first.coefficients.begin(),
                       result.coefficients.begin(),
                       std::minus<T>()); // l'ordre des arguments permet de ne parcourir que jusqu'au degré du plus petit polynôme.
    }
    result.adjust();
    result.dominant() = polynomial_ref_pair.second.dominant();
    return result;

    /*
     * The below alternative is slightly slower, (but much more readable...)
     */
//    Polynomial<T> result(q);
//    std::transform(result.coefficients.begin(), result.coefficients.end(), result.coefficients.begin(),
//                   std::negate<T>());
//    return p + result;
};

template<typename T>
Polynomial<T> operator*(const Polynomial<T> &p, const Polynomial<T> &q) {
    if (p.n < 0 || q.n < 0) {
        return Polynomial<T>();
    }
    int m = p.n + q.n;
    std::vector<T> coeffs(m + 1);
    for (int i = 0; i <= p.degree(); ++i) {
        for (int j = 0; j <= q.degree(); ++j) {
            coeffs[i + j] += p[i] * q[j];
        }
    }
    return Polynomial<T>(std::move(coeffs));
};

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

template<typename T>
Polynomial<T> operator/(const Polynomial<T> &p, const Polynomial<T> &q) {
    return euclid_div(p, q).first;
};

template<typename T>
Polynomial<T> operator%(const Polynomial<T> &p, const Polynomial<T> &q) {
    return euclid_div(p, q).second;
};

template<typename T>
Polynomial<T> operator+=(const Polynomial<T> &lhs, const Polynomial<T> &rhs) {
    return lhs + rhs;
};

template<typename T>
Polynomial<T> operator-=(const Polynomial<T> &lhs, const Polynomial<T> &rhs) {
    return lhs - rhs;
};

template<typename T>
Polynomial<T> operator*=(const Polynomial<T> &lhs, const Polynomial<T> &rhs) {
    return lhs * rhs;
};

template<typename T>
Polynomial<T> operator/=(const Polynomial<T> &lhs, const Polynomial<T> &rhs) {
    return lhs / rhs;
};

template<typename T>
Polynomial<T> operator%=(const Polynomial<T> &lhs, const Polynomial<T> &rhs) {
    return lhs % rhs;
}

template<typename T>
std::ostream &operator<<(std::ostream &out, const Polynomial<T> &p) {
    // +/- implementation ?
    if (p.n < 0) return out;
    for (int i = p.n; i > 0; --i) {
        if (!is_zero(p.coefficients[i])) {
            char var;

            if constexpr (std::is_same<T, std::complex<float>>::value || std::is_same<T, std::complex<double>>::value ||
                          std::is_same<T, std::complex<long double>>::value) { var = 'z'; }
            else { var = 'x'; }

            out << p.coefficients[i];

            if (i > 1) { out << "*" << var << "^" << i << " "; }
            else if (i == 1) { out << "*" << var; }
            out << " + ";
        }
    }
    if (!is_zero(p.coefficients[0]))
        out << p.coefficients[0];

    return out;
}

#endif //POLYNOMIALSHPC_POLYNOMIALS_HPP