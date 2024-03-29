// SparsePolynomials.hpp
// Created by Emmanuel Leclercq on 24/08/2023.
//

#ifndef POLYNOMIALSHPC_SPARSEPOLYNOMIALS_HPP
#define POLYNOMIALSHPC_SPARSEPOLYNOMIALS_HPP

#include <vector>
#include <list>
#include <map>
#include <unordered_map>
#include <queue>
#include <utility>
#include <algorithm>
#include "utils.hpp"
#include "Polynomials.hpp"

/*
 * Template declarations
 */
template<typename T>
class Monomial;

template<typename T1, typename T2>
auto operator<=>(const Monomial<T1> &u, const Monomial<T2> &v) {
    return u.degree() <=> v.degree();
}

template<typename T1, typename T2>
bool operator==(const Monomial<T1> &u, const Monomial<T2> &v) { return u.n == v.n && u.coeff() == v.coeff(); }

template<typename T1, typename T2>
bool operator==(const Monomial<T1> &u, const T2 &t) { return u.coeff() == t; }

template<typename T1, typename T2>
Monomial<decltype(T1() * T2())> operator*(const Monomial<T1> &lhs, const Monomial<T2> &rhs) {
    return Monomial<decltype(T1() * T2())>(lhs.coeff() * rhs.coeff(), lhs.degree() + rhs.degree());
}

template<typename T>
std::ostream &operator<<(std::ostream &, const Monomial<T> &);

template<typename T>
class SparsePolynomial;

template<typename T>
SparsePolynomial<T> derivative(SparsePolynomial<T> P, int k = 1);

template<typename T1, typename T2>
std::pair<SparsePolynomial<decltype(T1() * T2())>, SparsePolynomial<decltype(T1() * T2())>>
heap_div(const SparsePolynomial<T1> &, const SparsePolynomial<T2> &);

template<typename T>
std::ostream &operator<<(std::ostream &, const SparsePolynomial<T> &);

template<typename T1, typename T2>
SparsePolynomial<decltype(T1() * T2())> operator+(const SparsePolynomial<T1> &, const SparsePolynomial<T2> &);

template<typename T1, typename T2>
SparsePolynomial<decltype(T1() * T2())> operator-(const SparsePolynomial<T1> &, const SparsePolynomial<T2> &);

template<typename T1, typename T2>
SparsePolynomial<decltype(T1() * T2())> operator*(SparsePolynomial<T1> &, SparsePolynomial<T2> &);

template<typename T1, typename T2>
SparsePolynomial<decltype(T1() * T2())> operator/(const SparsePolynomial<T1> &, const SparsePolynomial<T2> &);

template<typename T1, typename T2>
SparsePolynomial<decltype(T1() * T2())> operator%(const SparsePolynomial<T1> &, const SparsePolynomial<T2> &);

template<typename T1, typename T2>
SparsePolynomial<decltype(T1() * T2())> operator+=(const SparsePolynomial<T2> &, const SparsePolynomial<T1> &);

template<typename T1, typename T2>
SparsePolynomial<decltype(T1() * T2())> operator-=(const SparsePolynomial<T1> &, const SparsePolynomial<T2> &);

template<typename T1, typename T2>
SparsePolynomial<decltype(T1() * T2())> operator*=(SparsePolynomial<T1> &, SparsePolynomial<T2> &);

template<typename T1, typename T2>
SparsePolynomial<decltype(T1() * T2())> operator/=(const SparsePolynomial<T1> &, const SparsePolynomial<T2> &);

template<typename T1, typename T2>
SparsePolynomial<decltype(T1() * T2())> operator%=(const SparsePolynomial<T1> &, const SparsePolynomial<T2> &);

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

    Monomial(const T &a, const int &m = 0) : coefficient(a), n(m) {}

    Monomial(const T &&a, const int &m = 0) : coefficient(a), n(m) {}

    /*
     * Ordering by degree
     */
    template<typename T1, typename T2>
    friend auto operator<=>(const Monomial<T1> &u, const Monomial<T2> &v);

    template<typename T1, typename T2>
    friend bool operator==(const Monomial<T1> &u, const Monomial<T2> &v);

    template<typename T1, typename T2>
    friend bool operator==(const Monomial<T1> &u, const T2 &t);

    friend std::ostream &operator
    <<<>(std::ostream &, const Monomial<T> &);

    template<typename T1, typename T2>
    friend Monomial<decltype(T1() * T2())> operator*(const Monomial<T1> &lhs, const Monomial<T2> &rhs);

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
    SparsePolynomial() : is_sorted(true) {}

    /*
     * Constructors assuming input is not 'sparse', which spears us
     * the need to know the degree beforehand
     */
    explicit SparsePolynomial(const std::vector<T> &coefficients);

    explicit SparsePolynomial(const std::vector<T> &&coefficients);

    explicit SparsePolynomial(const std::list<T> &coefficients);

    explicit SparsePolynomial(const std::list<T> &&coefficients);


    /* It makes more sense to use degrees as keys considering the
     * vectorial structure of polynomials spaces
     */
    explicit SparsePolynomial(const std::map<int, T, std::greater<int>> &coefficients);

    explicit SparsePolynomial(const std::map<int, T, std::greater<int>> &&coefficients);

    explicit SparsePolynomial(const std::map<int, T> &coefficients);

    explicit SparsePolynomial(const std::map<int, T> &&coefficients);

    explicit SparsePolynomial(const std::unordered_map<int, T> &coefficients);

    explicit SparsePolynomial(const std::unordered_map<int, T> &&coefficients);

    explicit SparsePolynomial(const std::vector<std::pair<int, T>> &coefficients);

    explicit SparsePolynomial(const std::vector<std::pair<int, T>> &&coefficients);

    explicit SparsePolynomial(const std::list<std::pair<int, T>> &coefficients);

    explicit SparsePolynomial(const std::list<std::pair<int, T>> &&coefficients);


    /*
     * Conversion from dense to sparse, with optional check for sparsity (more than half zeros)
     */
    explicit SparsePolynomial(const Polynomial<T> &P, bool check_sparsity = false);

    //  Default destructor will be used, standard data types
    ~SparsePolynomial() = default;

    [[nodiscard]] int degree() const { return n; };

    [[nodiscard]] Monomial<T> dominant() const { return *monomials.begin(); }

    [[nodiscard]] Monomial<T> &dominant() { return *monomials.begin(); }

    void reorder() {
        monomials.sort([](Monomial<T> &P, Monomial<T> &Q) { return Q < P; });
    }

    void adjust() {
        monomials.remove_if([](const Monomial<T> &M) { return M == 0; });
    }

    auto begin() const { return monomials.begin(); }

    auto end() const { return monomials.end(); }

    void add(const Monomial<T> &m) { if (m != 0) monomials.push_back(m); }

    template<typename U>
    T operator()(const U &) const;

    void derivative(int k = 1);

    friend SparsePolynomial<T> derivative(SparsePolynomial<T> P, int k);

    friend std::pair<SparsePolynomial<T>, SparsePolynomial<T>>
    euclid_div(const SparsePolynomial<T> &, const SparsePolynomial<T> &);

    template<typename T1, typename T2>
    friend std::pair<SparsePolynomial<decltype(T1() * T2())>, SparsePolynomial<decltype(T1() * T2())>>
    heap_div(const SparsePolynomial<T1> &, const SparsePolynomial<T2> &);

    friend std::ostream &operator
    <<<>(std::ostream &, const SparsePolynomial<T> &);

    template<typename T1, typename T2>
    friend SparsePolynomial<decltype(T1() * T2())>
    operator+(const SparsePolynomial<T1> &, const SparsePolynomial<T2> &);

    template<typename T1, typename T2>
    friend SparsePolynomial<decltype(T1() * T2())>
    operator-(const SparsePolynomial<T1> &, const SparsePolynomial<T2> &);

    template<typename T1, typename T2>
    friend SparsePolynomial<decltype(T1() * T2())>
    operator*(SparsePolynomial<T1> &, SparsePolynomial<T2> &);

    template<typename T1, typename T2>
    friend SparsePolynomial<decltype(T1() * T2())>
    operator/(const SparsePolynomial<T1> &, const SparsePolynomial<T2> &);

    template<typename T1, typename T2>
    friend SparsePolynomial<decltype(T1() * T2())>
    operator%(const SparsePolynomial<T1> &, const SparsePolynomial<T2> &);

    template<typename T1, typename T2>
    friend SparsePolynomial<decltype(T1() * T2())>
    operator+=(const SparsePolynomial<T1> &, const SparsePolynomial<T2> &);

    template<typename T1, typename T2>
    friend SparsePolynomial<decltype(T1() * T2())>
    operator-=(const SparsePolynomial<T1> &, const SparsePolynomial<T2> &);

    template<typename T1, typename T2>
    friend SparsePolynomial<decltype(T1() * T2())>
    operator*=(SparsePolynomial<T1> &, SparsePolynomial<T2> &);

    template<typename T1, typename T2>
    friend SparsePolynomial<decltype(T1() * T2())>
    operator/=(const SparsePolynomial<T1> &, const SparsePolynomial<T2> &);

    template<typename T1, typename T2>
    friend SparsePolynomial<decltype(T1() * T2())>
    operator%=(const SparsePolynomial<T1> &, const SparsePolynomial<T2> &);
};

/*
 * ! Zeros are stored too, maybe we can strip them immediately here,
 * without additional complexity?
*/
template<typename T>
SparsePolynomial<T>::SparsePolynomial(const std::vector<T> &coefficients) {
    int index = coefficients.size();
    n = index;
    std::ranges::for_each(coefficients | std::views::reverse, [&index, this](const T &p) {
        index--;
        if (p != 0) { monomials.push_back(Monomial<T>(p, index)); }
    });
}

template<typename T>
SparsePolynomial<T>::SparsePolynomial(const std::vector<T> &&coefficients) {
    int index = coefficients.size();
    n = index;
    std::ranges::for_each(coefficients | std::views::reverse, [&index, this](const T &p) {
        index--;
        if (p != 0) { monomials.push_back(Monomial<T>(std::move(p), index)); }
    });
}

template<typename T>
SparsePolynomial<T>::SparsePolynomial(const std::list<T> &coefficients) {
    int index = coefficients.size();
    n = index;
    std::ranges::for_each(coefficients | std::views::reverse, [&index, this](const T &p) {
        index--;
        if (p != 0) { monomials.push_back(Monomial<T>(p, index)); }
    });
}

template<typename T>
SparsePolynomial<T>::SparsePolynomial(const std::list<T> &&coefficients) {
    int index = coefficients.size();
    n = index;
    std::ranges::for_each(coefficients | std::views::reverse, [&index, this](const T &p) {
        index--;
        if (p != 0) { monomials.push_back(Monomial<T>(p, index)); }
    });
}

template<typename T>
SparsePolynomial<T>::SparsePolynomial(const std::map<int, T, std::greater<int>> &coefficients) {
    std::ranges::transform(coefficients, std::back_inserter(this->monomials),
                           [this](const std::pair<int, T> &p) {
                               n = std::max(n, p.first);
                               return Monomial<T>(p.second, p.first);
                           });
}

template<typename T>
SparsePolynomial<T>::SparsePolynomial(const std::map<int, T, std::greater<int>> &&coefficients) {
    std::ranges::transform(coefficients, std::back_inserter(this->monomials),
                           [this](std::pair<int, T> &p) {
                               n = std::max(n, p.first);
                               return Monomial<T>(std::move(p.second), p.first);
                           });
}

template<typename T>
SparsePolynomial<T>::SparsePolynomial(const std::map<int, T> &coefficients) {
    std::ranges::transform(coefficients, std::front_inserter(this->monomials),
                           [this](const std::pair<int, T> &p) {
                               n = std::max(n, p.first);
                               return Monomial<T>(p.second, p.first);
                           });
}

template<typename T>
SparsePolynomial<T>::SparsePolynomial(const std::map<int, T> &&coefficients) {
    std::ranges::transform(coefficients.begin(), coefficients.end(), std::front_inserter(this->monomials),
                           [this](std::pair<int, T> &p) {
                               n = std::max(n, p.first);
                               return Monomial<T>(std::move(p.second), p.first);
                           });
}

template<typename T>
SparsePolynomial<T>::SparsePolynomial(const std::unordered_map<int, T> &&coefficients) {
    std::ranges::transform(coefficients, std::back_inserter(this->monomials),
                           [this](std::pair<int, T> &p) {
                               n = std::max(n, p.first);
                               return Monomial<T>(std::move(p.second), p.first);
                           });
    this->reorder();
}

template<typename T>
SparsePolynomial<T>::SparsePolynomial(const std::unordered_map<int, T> &coefficients) {
    /*
     * We could the unordered_map elements to a heap after transformation and then feed
     * the monomials list so that it comes out sorted, instead of sorting separately.
     * Might be worth testing which is faster.
     */
    std::ranges::transform(coefficients, std::back_inserter(this->monomials), [this](const std::pair<int, T> p) {
        n = std::max(n, p.first);
        return Monomial<T>(p.second, p.first);
    });
    this->reorder();
}


template<typename T>
SparsePolynomial<T>::SparsePolynomial(const std::vector<std::pair<int, T>> &coefficients) {
    /*
     * The input vector is assumed ordered by degree
     */
    std::ranges::transform(coefficients, std::back_inserter(this->monomials),
                           [this](const std::pair<int, T> &p) {
                               n = std::max(n, p.first);
                               return Monomial<T>(p.second, p.first);
                           });
}

template<typename T>
SparsePolynomial<T>::SparsePolynomial(const std::vector<std::pair<int, T>> &&coefficients) {
    /*
     * The input vector is assumed ordered by degree
     */
    std::ranges::transform(coefficients, std::back_inserter(this->monomials),
                           [this](std::pair<int, T> &p) {
                               n = std::max(n, p.first);
                               return Monomial<T>(std::move(p.second), p.first);
                           });
}

template<typename T>
SparsePolynomial<T>::SparsePolynomial(const std::list<std::pair<int, T>> &coefficients) {
    /*
     * The input list is assumed ordered by degree
     */
    std::ranges::transform(coefficients, std::back_inserter(this->monomials),
                           [this](const std::pair<int, T> &p) {
                               n = std::max(n, p.first);
                               return Monomial<T>(p.second, p.first);
                           });
}

template<typename T>
SparsePolynomial<T>::SparsePolynomial(const std::list<std::pair<int, T>> &&coefficients) {
    /*
     * The input list is assumed ordered by degree
     */
    std::ranges::transform(coefficients, std::back_inserter(this->monomials), [this](std::pair<int, T> &p) {
        n = std::max(n, p.first);
        return Monomial<T>(std::move(p.second), p.first);
    });
}

template<typename T>
SparsePolynomial<T>::SparsePolynomial(const Polynomial<T> &P, bool check_sparsity) {
    if (check_sparsity && !P.is_sparse()) {
        throw std::invalid_argument(
                "The polynomial is not sparse enough, to force conversion set second parameter to false");
    }

    n = P.degree();
    for (int index = n; index > -1; index--) {
        if (!is_zero(P[index])) {
            monomials.emplace_back(P[index], index);
        }
    }
}

template<typename T>
void SparsePolynomial<T>::derivative(int k) {
    if (k == 0) { return; }
    if (degree() < k) {
        *this = SparsePolynomial<T>();
        return;
    }
    if (k < 0) {
        for (auto it = monomials.begin(); it != monomials.end(); it++) {
            it->coeff() /= rangeProduct<T>(it->degree() + 1, it->degree() - k);
            it->degree() -= k;
        }
        n -= k;
    } else {
        while (monomials.back().degree() < k) {
            monomials.pop_back();
        }
        for (auto it = monomials.begin(); it != monomials.end(); it++) {
            it->coeff() *= rangeProduct<T>(it->degree() - k + 1, it->degree());
            it->degree() -= k;
        }
        n -= k;
    }
}

template<typename T>
SparsePolynomial<T> derivative(SparsePolynomial<T> P, int k) {
    P.derivative(k);
    return P;
}

template<typename T>
template<typename U>
T SparsePolynomial<T>::operator()(const U &x) const {
    T result = 0;
    for (const auto &m: monomials) {
        result += m.coeff() * std::pow(x, m.degree());
    }
    return result;
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

template<typename T1, typename T2>
SparsePolynomial<decltype(T1() * T2())> operator+(const SparsePolynomial<T1> &lhs, const SparsePolynomial<T2> &rhs) {
    SparsePolynomial<decltype(T1() * T2())> result;
    auto it1 = lhs.begin();
    auto it2 = rhs.begin();
    while (it1 != lhs.end() && it2 != rhs.end()) {
        if (it1->degree() > it2->degree()) {
            result.monomials.push_back(*it1);
            ++it1;
        } else if (it1->degree() < it2->degree()) {
            result.monomials.push_back(*it2);
            ++it2;
        } else {
            if (decltype(T1() * T2()) new_coefficient = it1->coeff() + it2->coeff();
                    new_coefficient != decltype(T1() * T2())(0)) {
                result.monomials.push_back(Monomial<decltype(T1() * T2())>(it1->degree(), new_coefficient));
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
}

template<typename T1, typename T2>
SparsePolynomial<decltype(T1() * T2())> operator-(const SparsePolynomial<T1> &lhs, const SparsePolynomial<T2> &rhs) {
    SparsePolynomial<decltype(T1() * T2())> result(rhs);
    std::ranges::transform(result.monomials, begin(result.monomials),
                           [](const Monomial<decltype(T1() * T2())> &m) { return negate(m); });
    return lhs + result;
}

template<typename T1, typename T2>
SparsePolynomial<decltype(T1() * T2())> operator*(SparsePolynomial<T1> &P, SparsePolynomial<T2> &Q) {
    std::priority_queue<Monomial<decltype(T1() * T2())>> heap;
    SparsePolynomial<decltype(T1() * T2())> result;
    for (const auto &monoL: P) {
        for (const auto &monoR: Q) {
            heap.push(monoL * monoR);
        }
    }
    while (!heap.empty()) {
        auto mono = heap.top();
        heap.pop();
        while (!heap.empty() && heap.top().degree() == mono.degree()) {
            mono.coeff() += heap.top().coeff();
            heap.pop();
        }
        result.add(mono);
    }

    return result;

}

template<typename T>
std::pair<SparsePolynomial<T>, SparsePolynomial<T>>
euclid_div(const SparsePolynomial<T> &P, const SparsePolynomial<T> &Q) {
    SparsePolynomial<T> quotient;
    auto remainder = P;
    while (!remainder.monomials.empty() && remainder.monomials.back().degree() >= Q.monomials.back().degree()) {

    }

    return std::pair<SparsePolynomial<T>, SparsePolynomial<T>>(quotient, remainder);
}

template<typename T1, typename T2>
std::pair<SparsePolynomial<decltype(T1() * T2())>, SparsePolynomial<decltype(T1() * T2())>>
heap_div(const SparsePolynomial<T1> &dividend, const SparsePolynomial<T2> &divisor) {
    SparsePolynomial<decltype(T1() * T2())> quotient;
    SparsePolynomial<decltype(T1() * T2())> remainder = dividend;
    std::priority_queue<Monomial<decltype(T1() * T2())>> heap;
    for (const auto &mono: dividend) {
        heap.push(mono);
    }
    while (!heap.empty() && heap.top().degree() >= divisor.degree()) {
        auto LeadingTermDividend = heap.top();
        Monomial<decltype(T1() * T2())> QuotientTerm{LeadingTermDividend.coeff() / divisor.dominant().coeff(),
                                                     LeadingTermDividend.degree() - divisor.dominant().degree()};
        quotient.add(QuotientTerm);

    }
    return {quotient, remainder};
}

template<typename T1, typename T2>
SparsePolynomial<decltype(T1() * T2())> operator/(const SparsePolynomial<T1> &P, const SparsePolynomial<T2> &Q) {
    SparsePolynomial<decltype(T1() * T2())> quotient;
    auto remainder = P;

    while (!remainder.monomials.empty() && remainder.monomials.back().degree() >= Q.monomials.back().degree()) {
        auto leadCoeff = remainder.monomials.back().coeff() / Q.monomials.back().coeff();
        int leadDegree = remainder.monomials.back().degree() - Q.monomials.back().degree();
        quotient.monomials.push_back(Monomial<decltype(T1() * T2())>{leadCoeff, leadDegree});

        SparsePolynomial<decltype(T1() * T2())> product;
        for (const auto &mono: Q) {
            product.monomials.push_back(
                    Monomial<decltype(T1() * T2())>{mono.coeff() * leadCoeff, mono.degree() + leadDegree});
        }

        remainder = remainder - product;
    }
    return quotient;
}

template<typename T1, typename T2>
SparsePolynomial<decltype(T1() * T2())> operator%(const SparsePolynomial<T1> &, const SparsePolynomial<T2> &) {

}

template<typename T1, typename T2>
SparsePolynomial<decltype(T1() * T2())> operator+=(const SparsePolynomial<T1> &A, const SparsePolynomial<T2> &B) {
    return A + B;
}

template<typename T1, typename T2>
SparsePolynomial<decltype(T1() * T2())> operator-=(const SparsePolynomial<T1> &A, const SparsePolynomial<T2> &B) {
    return A - B;
}

template<typename T1, typename T2>
SparsePolynomial<decltype(T1() * T2())> operator*=(SparsePolynomial<T1> &A, SparsePolynomial<T2> &B) {
    return A * B;
}

template<typename T1, typename T2>
SparsePolynomial<decltype(T1() * T2())>
operator/=(const SparsePolynomial<T1> &, const SparsePolynomial<T2> &) {/* TODO */}

template<typename T1, typename T2>
SparsePolynomial<decltype(T1() * T2())>
operator%=(const SparsePolynomial<T1> &, const SparsePolynomial<T2> &) {/* TODO */}

#endif //POLYNOMIALSHPC_SPARSEPOLYNOMIALS_HPP
