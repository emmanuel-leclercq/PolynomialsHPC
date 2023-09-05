//
// Created by Emmanuel Leclercq on 24/08/2023.
//

#ifndef POLYNOMIALSHPC_SPARSEPOLYNOMIALS_HPP
#define POLYNOMIALSHPC_SPARSEPOLYNOMIALS_HPP

#include <vector>
#include <list>
#include <map>
#include <unordered_map>
#include <algorithm>
#include "utils.hpp"
#include "Polynomials.hpp"

/*
 * Template declarations
 */
template<typename T>
class Monomial;

template<typename T, typename Y>
bool operator<(const Monomial<T> &u, const Monomial<Y> &v) { return u.n < v.n; }

template<typename T, typename Y>
bool operator==(const Monomial<T> &u, const Monomial<Y> &v) { return u.n == v.n && u.coeff() == v.coeff(); }


template<typename T, typename Y>
bool operator==(const Monomial<T> &u, const Y &t) { return u.coeff() == t; }

template<typename T>
std::ostream &operator<<(std::ostream &, const Monomial<T> &);

template<typename T>
class SparsePolynomial;

template<typename T>
std::ostream &operator<<(std::ostream &, const SparsePolynomial<T> &);

template<typename T>
SparsePolynomial<T> operator+(const SparsePolynomial<T> &, const SparsePolynomial<T> &);

template<typename T>
SparsePolynomial<T> operator-(const SparsePolynomial<T> &, const SparsePolynomial<T> &);

template<typename T>
SparsePolynomial<T> operator*(const SparsePolynomial<T> &, const SparsePolynomial<T> &);

template<typename T>
SparsePolynomial<T> operator/(const SparsePolynomial<T> &, const SparsePolynomial<T> &);

template<typename T>
SparsePolynomial<T> operator%(const SparsePolynomial<T> &, const SparsePolynomial<T> &);

template<typename T>
SparsePolynomial<T> operator+=(const SparsePolynomial<T> &, const SparsePolynomial<T> &);

template<typename T>
SparsePolynomial<T> operator-=(const SparsePolynomial<T> &, const SparsePolynomial<T> &);

template<typename T>
SparsePolynomial<T> operator*=(const SparsePolynomial<T> &, const SparsePolynomial<T> &);

template<typename T>
SparsePolynomial<T> operator/=(const SparsePolynomial<T> &, const SparsePolynomial<T> &);

template<typename T>
SparsePolynomial<T> operator%=(const SparsePolynomial<T> &, const SparsePolynomial<T> &);


/*
 * Definition of the parent class monomial and all its methods
 */
template<typename T>
class Monomial {
private:
    int n;
    T coefficient;

public:
    Monomial() : coefficient(0), n(-1) {}

    explicit Monomial(const T &a, const int &m = 0) : coefficient(a), n(m) {}

    explicit Monomial(const T &&a, const int &m = 0) : coefficient(a), n(m) {}

    /*
     * Ordering by degree
     */
    friend bool operator< /* */ <>(const Monomial<T> &u, const Monomial<T> &v);

    friend bool operator==<>(const Monomial<T> &u, const Monomial<T> &v);

    friend bool operator==<>(const Monomial<T> &u, const T &t);

    friend std::ostream &operator
    <<<>(std::ostream &, const Monomial<T> &);

    T coeff() const { return coefficient; }

    [[nodiscard]] int degree() const { return n; }

    T &coeff() { return coefficient; }

    int &degree() { return n; }
};

template<typename T>
std::ostream &operator<<(std::ostream &out, const Monomial<T> &p) {

    if (p.n == 0) [[unlikely]] {
        out << p.coefficient;
    } else if (is_zero(p.coefficient)) {
        out << "0";
    } else if (is_one(p.coefficient)) {
        if (p.n != 1)[[likely]] { out << "x^" << p.n; }
        else { out << "x"; }
    } else {
        if (p.n != 1)[[likely]] { out << p.coefficient << "x^" << p.n; }
        else { out << p.coefficient << "x"; }
    }
    return out;
}

template<typename T>
Monomial<T> negate(Monomial<T> m) { return Monomial<T>(-m.coeff(), m.degree()); }

/*
 * Definition of the Sparse Polynomial class and all of its methods
 */

template<typename T>
class SparsePolynomial {
private:
    std::list<Monomial<T>> monomials;
    int n = -1;
    bool is_sorted = false;

public:
    /*
     * Providing a handful of constructors for convenience and flexibility,
     * considering certain data types make more sense for a sparse implementation,
     * and to reduce the need for preprocessing existing variables
     */
    SparsePolynomial() : monomials(), n(-1), is_sorted(true) {}

    /*
     * Constructors assuming input is not 'sparse', which spears us
     * the need to know the degree beforehand
     */
    explicit SparsePolynomial(const std::vector<T> &coefficients, bool sorted);

    explicit SparsePolynomial(const std::vector<T> &&coefficients, bool sorted);

    explicit SparsePolynomial(const std::list<T> &coefficients, bool sorted);

    explicit SparsePolynomial(const std::list<T> &&coefficients, bool sorted);


    /* It makes more sense to use degrees as keys considering the
     * vectorial structure of polynomials spaces
     */
    explicit SparsePolynomial(const std::map<int, T, std::greater<int>> &coefficients, bool sorted);

    explicit SparsePolynomial(const std::map<int, T, std::greater<int>> &&coefficients, bool sorted);

    explicit SparsePolynomial(const std::map<int, T> &coefficients, bool sorted);

    explicit SparsePolynomial(const std::map<int, T> &&coefficients, bool sorted);

    explicit SparsePolynomial(const std::unordered_map<int, T> &coefficients, bool sorted);

    explicit SparsePolynomial(const std::unordered_map<int, T> &&coefficients, bool sorted);

    explicit SparsePolynomial(const std::vector<std::pair<int, T>> &coefficients, bool sorted);

    explicit SparsePolynomial(const std::vector<std::pair<int, T>> &&coefficients, bool sorted);

    explicit SparsePolynomial(const std::list<std::pair<int, T>> &coefficients, bool sorted);

    explicit SparsePolynomial(const std::list<std::pair<int, T>> &&coefficients, bool sorted);


    /*
     * Conversion from dense to sparse, with optional check for sparsity (more than half zeros)
     */
    explicit SparsePolynomial(const Polynomial<T> &P, bool check_sparsity = false);

    //  Default destructor will be used, standard data types
    ~SparsePolynomial() = default;

    [[nodiscard]] int degree() const { return n; }

    [[nodiscard]] Monomial<T> dominant() const { return *monomials.begin(); }

    [[nodiscard]] Monomial<T> &dominant() { return *monomials.begin(); }

    void reorder() {
        monomials.sort([](Monomial<T> &P, Monomial<T> &Q) { return Q < P; });
        is_sorted = true;
    }

    void adjust() {
        monomials.remove_if([](const Monomial<T> &M) { return M == 0; });
    }

    auto begin() const { return monomials.begin(); }

    auto end() const { return monomials.end(); }

    friend std::ostream &operator
    <<<>(std::ostream &, const SparsePolynomial<T> &);

    friend SparsePolynomial<T> operator+<T>(const SparsePolynomial<T> &, const SparsePolynomial<T> &);

    friend SparsePolynomial<T> operator-<T>(const SparsePolynomial<T> &, const SparsePolynomial<T> &);

    friend SparsePolynomial<T> operator*<T>(const SparsePolynomial<T> &, const SparsePolynomial<T> &);

    friend SparsePolynomial<T> operator/<T>(const SparsePolynomial<T> &, const SparsePolynomial<T> &);

    friend SparsePolynomial<T> operator%<T>(const SparsePolynomial<T> &, const SparsePolynomial<T> &);

    friend SparsePolynomial<T> operator+=<T>(const SparsePolynomial<T> &, const SparsePolynomial<T> &);

    friend SparsePolynomial<T> operator-=<T>(const SparsePolynomial<T> &, const SparsePolynomial<T> &);

    friend SparsePolynomial<T> operator*=<T>(const SparsePolynomial<T> &, const SparsePolynomial<T> &);

    friend SparsePolynomial<T> operator/=<T>(const SparsePolynomial<T> &, const SparsePolynomial<T> &);

    friend SparsePolynomial<T> operator%=<T>(const SparsePolynomial<T> &, const SparsePolynomial<T> &);
};

/*
 * ! Zeros are stored too, maybe we can strip them immediately here,
 * without additional complexity?
*/
template<typename T>
SparsePolynomial<T>::SparsePolynomial(const std::vector<T> &coefficients, bool sorted) {
    int index = coefficients.size();
    n = index;
    std::transform(coefficients.begin(), coefficients.end(), std::back_inserter(this->monomials),
                   [&index](const T &p) { return Monomial<T>(p, index--); });
}

template<typename T>
SparsePolynomial<T>::SparsePolynomial(const std::vector<T> &&coefficients, bool sorted) {
    int index = coefficients.size();
    n = index;
    std::transform(coefficients.begin(), coefficients.end(), std::back_inserter(this->monomials),
                   [&index](T &p) { return Monomial<T>(std::move(p), index--); });
}

template<typename T>
SparsePolynomial<T>::SparsePolynomial(const std::list<T> &coefficients, bool sorted) {
    int index = coefficients.size();
    n = index;
    std::transform(coefficients.begin(), coefficients.end(), std::back_inserter(this->monomials),
                   [&index](const T &p) { return Monomial<T>(p, index--); });
}

template<typename T>
SparsePolynomial<T>::SparsePolynomial(const std::list<T> &&coefficients, bool sorted) {
    int index = coefficients.size();
    n = index;
    std::transform(coefficients.begin(), coefficients.end(), std::back_inserter(this->monomials),
                   [&index](T &p) { return Monomial<T>(std::move(p), index--); });
}

template<typename T>
SparsePolynomial<T>::SparsePolynomial(const std::map<int, T, std::greater<int>> &coefficients, bool sorted) {
    std::transform(coefficients.begin(), coefficients.end(), std::back_inserter(this->monomials),
                   [this](const std::pair<int, T> &p) {
                       n = std::max(n, p.first);
                       return Monomial<T>(p.second, p.first);
                   });
    this->is_sorted = true;
}

template<typename T>
SparsePolynomial<T>::SparsePolynomial(const std::map<int, T, std::greater<int>> &&coefficients, bool sorted) {
    std::transform(coefficients.begin(), coefficients.end(), std::back_inserter(this->monomials),
                   [this](std::pair<int, T> &p) {
                       n = std::max(n, p.first);
                       return Monomial<T>(std::move(p.second), p.first);
                   });
    this->is_sorted = true;
    delete coefficients;
}

template<typename T>
SparsePolynomial<T>::SparsePolynomial(const std::map<int, T> &coefficients, bool sorted) {
    std::transform(coefficients.begin(), coefficients.end(), std::back_inserter(this->monomials),
                   [this](const std::pair<int, T> &p) {
                       n = std::max(n, p.first);
                       return Monomial<T>(p.second, p.first);
                   });
}

template<typename T>
SparsePolynomial<T>::SparsePolynomial(const std::map<int, T> &&coefficients, bool sorted) {
    std::transform(coefficients.begin(), coefficients.end(), std::back_inserter(this->monomials),
                   [this](std::pair<int, T> &p) {
                       n = std::max(n, p.first);
                       return Monomial<T>(std::move(p.second), p.first);
                   });

    delete coefficients;
}

template<typename T>
SparsePolynomial<T>::SparsePolynomial(const std::unordered_map<int, T> &&coefficients, bool sorted) {
    std::transform(coefficients.begin(), coefficients.end(), std::back_inserter(this->monomials),
                   [this](std::pair<int, T> &p) {
                       n = std::max(n, p.first);
                       return Monomial<T>(std::move(p.second), p.first);
                   });
    delete coefficients;
}

template<typename T>
SparsePolynomial<T>::SparsePolynomial(const std::unordered_map<int, T> &coefficients, bool sorted) {
    std::transform(coefficients.begin(), coefficients.end(), std::back_inserter(this->monomials),
                   [this](const std::pair<int, T> &p) {
                       n = std::max(n, p.first);
                       return Monomial<T>(p.second, p.first);
                   });
}


template<typename T>
SparsePolynomial<T>::SparsePolynomial(const std::vector<std::pair<int, T>> &coefficients, bool sorted) {
    std::transform(coefficients.begin(), coefficients.end(), std::back_inserter(this->monomials),
                   [this](const std::pair<int, T> &p) {
                       n = std::max(n, p.first);
                       return Monomial<T>(p.second, p.first);
                   });
}

template<typename T>
SparsePolynomial<T>::SparsePolynomial(const std::vector<std::pair<int, T>> &&coefficients, bool sorted) {
    std::transform(coefficients.begin(), coefficients.end(), std::back_inserter(this->monomials),
                   [this](std::pair<int, T> &p) {
                       n = std::max(n, p.first);
                       return Monomial<T>(std::move(p.second), p.first);
                   });
    delete coefficients;
}

template<typename T>
SparsePolynomial<T>::SparsePolynomial(const std::list<std::pair<int, T>> &coefficients, bool sorted) {
    std::transform(coefficients.begin(), coefficients.end(), std::back_inserter(this->monomials),
                   [this](const std::pair<int, T> &p) {
                       n = std::max(n, p.first);
                       return Monomial<T>(p.second, p.first);
                   });
}

template<typename T>
SparsePolynomial<T>::SparsePolynomial(const std::list<std::pair<int, T>> &&coefficients, bool sorted) {
    std::transform(coefficients.begin(), coefficients.end(), std::back_inserter(this->monomials),
                   [this](std::pair<int, T> &p) {
                       n = std::max(n, p.first);
                       return Monomial<T>(std::move(p.second), p.first);
                   });
    delete coefficients;
}

template<typename T>
SparsePolynomial<T>::SparsePolynomial(const Polynomial<T> &P, bool check_sparsity) {
    if (check_sparsity) {
        if (!P.is_sparse()) {
            throw std::invalid_argument("The polynomial is not sparse enough");
        }
    }
    n = P.degree();
    for (int index = n; index > -1; index--) {
        if (!is_zero(P[index])) {
            monomials.emplace_back(P[index], index);
        }
    }
}

template<typename T>
std::ostream &operator<<(std::ostream &out, const SparsePolynomial<T> &P) {
    for (auto it = P.monomials.begin(); it != P.monomials.end(); ++it) {
        out << *it;
        if (it != P.monomials.end() && std::next(it) != P.monomials.end()) {
            if (std::next(it)->coeff() >= 0) { out << " +"; }
            else { out << " "; }
        }
    }
    return out;
}

template<typename T>
SparsePolynomial<T> operator+(const SparsePolynomial<T> &lhs, const SparsePolynomial<T> &rhs) {
    if (rhs.is_sorted && lhs.is_sorted) {
        SparsePolynomial<T> result;
        auto it1 = lhs.begin();
        auto it2 = rhs.begin();
        while (it1 != lhs.end() && it2 != rhs.end()) {
            if (it1->degree() < it2->degree()) {
                result.monomials.push_back(*it1);
                ++it1;
            } else if (it1->degree() > it2->degree()) {
                result.monomials.push_back(*it2);
                ++it2;
            } else {
                T new_coefficient = it1->coeff() + it2->coeff();

                if (new_coefficient != T(0)) {
                    result.monomials.push_back(Monomial<T>(it1->degree(), new_coefficient));
                }

                ++it1;
                ++it2;
            }
        }

        // Handle the remaining monomials in either polynomial
        while (it1 != lhs.end()) {
            result.monomials.push_back(*it1);
            ++it1;
        }
        while (it2 != rhs.end()) {
            result.monomials.push_back(*it2);
            ++it2;
        }
        return result;
    } else {
        std::map<int, T, std::greater<int>> result;
        for (auto it = lhs.begin(); it != lhs.end(); it++) {
            result[it->degree()] += it->coeff();
        }
        for (auto it = rhs.begin(); it != rhs.end(); it++) {
            result[it->degree()] += it->coeff();
        }
        SparsePolynomial<T> ans(result, true);
        return ans;
    }
}

template<typename T>
SparsePolynomial<T> operator-(const SparsePolynomial<T> &lhs, const SparsePolynomial<T> &rhs) {
    SparsePolynomial<T> result(rhs);
    std::transform(result.monomials.begin(), result.monomials.end(), result.monomials.begin(),
                   [](const Monomial<T> &m) { return negate(m); });
    return lhs + result;
}

template<typename T>
SparsePolynomial<T> operator*(const SparsePolynomial<T> &P, const SparsePolynomial<T> &Q) {
    if (P.is_sorted && Q.is_sorted) {

    } else {

    }
}

template<typename T>
SparsePolynomial<T> operator/(const SparsePolynomial<T> &, const SparsePolynomial<T> &) {}

template<typename T>
SparsePolynomial<T> operator%(const SparsePolynomial<T> &, const SparsePolynomial<T> &) {}

template<typename T>
SparsePolynomial<T> operator+=(const SparsePolynomial<T> &, const SparsePolynomial<T> &) {}

template<typename T>
SparsePolynomial<T> operator-=(const SparsePolynomial<T> &, const SparsePolynomial<T> &) {}

template<typename T>
SparsePolynomial<T> operator*=(const SparsePolynomial<T> &, const SparsePolynomial<T> &) {}

template<typename T>
SparsePolynomial<T> operator/=(const SparsePolynomial<T> &, const SparsePolynomial<T> &) {}

template<typename T>
SparsePolynomial<T> operator%=(const SparsePolynomial<T> &, const SparsePolynomial<T> &) {}

#endif //POLYNOMIALSHPC_SPARSEPOLYNOMIALS_HPP
