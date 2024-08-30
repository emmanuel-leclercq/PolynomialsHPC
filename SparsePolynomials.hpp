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
#include "DensePolynomials.hpp"

namespace Polynomial {
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
    class Sparse;

    template<typename T>
    Sparse<T> derivative(Sparse<T> P, int k = 1);

    template<typename T1, typename T2>
    std::pair<Sparse<decltype(T1() * T2())>, Sparse<decltype(T1() * T2())>>
    heap_div(const Sparse<T1> &, const Sparse<T2> &);

    template<typename T>
    std::ostream &operator<<(std::ostream &, const Sparse<T> &);

    template<typename T1, typename T2>
    bool operator==(const Sparse<T1> &, const Sparse<T2> &);

    template<typename T1, typename T2>
    Sparse<decltype(T1() * T2())> operator+(const Sparse<T1> &, const Sparse<T2> &);

    template<typename T1, typename T2>
    Sparse<decltype(T1() * T2())> operator-(const Sparse<T1> &, const Sparse<T2> &);

    template<typename T1, typename T2>
    Sparse<decltype(T1() * T2())> operator*(Sparse<T1> &, Sparse<T2> &);

    template<typename T1, typename T2>
    Sparse<decltype(T1() * T2())> operator/(const Sparse<T1> &, const Sparse<T2> &);

    template<typename T1, typename T2>
    Sparse<decltype(T1() * T2())> operator%(const Sparse<T1> &, const Sparse<T2> &);

    template<typename T1, typename T2>
    Sparse<decltype(T1() * T2())> operator+=(const Sparse<T2> &, const Sparse<T1> &);

    template<typename T1, typename T2>
    Sparse<decltype(T1() * T2())> operator-=(const Sparse<T1> &, const Sparse<T2> &);

    template<typename T1, typename T2>
    Sparse<decltype(T1() * T2())> operator*=(Sparse<T1> &, Sparse<T2> &);

    template<typename T1, typename T2>
    Sparse<decltype(T1() * T2())> operator/=(const Sparse<T1> &, const Sparse<T2> &);

    template<typename T1, typename T2>
    Sparse<decltype(T1() * T2())> operator%=(const Sparse<T1> &, const Sparse<T2> &);

/*
 * Definition of the parent class monomial and all its methods
 */
    template<typename T>
    class Monomial {
    private:
        int n;
        T coefficient;

    public:
        Monomial() : n(-1), coefficient(0) {}

        explicit Monomial(const T &a, const int &m = 0) : n(m), coefficient(a) { if (a == 0) { n = -1; }}

        explicit Monomial(const T &&a, const int &m = 0) : n(m), coefficient(a) { if (a == 0) { n = -1; }}

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
    class Sparse {
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
        Sparse() : is_sorted(true) {}

        /*
         * Constructors assuming input is not 'sparse', which spears us
         * the need to know the degree beforehand
         */
        explicit Sparse(const std::vector<T> &coefficients);

        explicit Sparse(const std::vector<T> &&coefficients);

        explicit Sparse(const std::list<T> &coefficients);

        explicit Sparse(const std::list<T> &&coefficients);


        /* It makes more sense to use degrees as keys considering the
         * vectorial structure of polynomials spaces
         */
        explicit Sparse(const std::map<int, T, std::greater<>> &coefficients);

        explicit Sparse(const std::map<int, T, std::greater<>> &&coefficients);

        explicit Sparse(const std::map<int, T> &coefficients);

        explicit Sparse(const std::map<int, T> &&coefficients);

        explicit Sparse(const std::unordered_map<int, T> &coefficients);

        explicit Sparse(const std::unordered_map<int, T> &&coefficients);

        explicit Sparse(const std::vector<std::pair<int, T>> &coefficients);

        explicit Sparse(const std::vector<std::pair<int, T>> &&coefficients);

        explicit Sparse(const std::list<std::pair<int, T>> &coefficients);

        explicit Sparse(const std::list<std::pair<int, T>> &&coefficients);


        /*
         * Conversion from dense to sparse, with optional check for sparsity (more than half zeros)
         */
        explicit Sparse(const Polynomial::Dense<T> &P, bool check_sparsity = false);

        //  Default destructor will be used, standard data types
        ~Sparse() = default;

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

        /*
         * Accessing the Monomials by degree
         */
        template<std::integral n>
        Monomial<T> operator[](n i) const {
            return *std::ranges::find_if(monomials, [i](const Monomial<T> &m) { return m.degree() == i; });
        }

        void derivative(int k = 1);

        template<typename U>
        friend Sparse<U> derivative(Sparse<U> P, int k);

        friend std::pair<Sparse<T>, Sparse<T>>
        euclid_div(const Sparse<T> &, const Sparse<T> &);

        template<typename T1, typename T2>
        friend std::pair<Sparse<decltype(T1() * T2())>, Sparse<decltype(T1() * T2())>>
        heap_div(const Sparse<T1> &, const Sparse<T2> &);

        friend std::ostream &operator
        <<<>(std::ostream &, const Sparse<T> &);

        template<typename T1, typename T2>
        friend bool operator==(const Sparse<T1> &, const Sparse<T2> &);

        bool operator!=(const Sparse &rhs) const {
            return !(rhs == *this);
        }

        template<typename T1, typename T2>
        friend Sparse<decltype(T1() * T2())>
        operator+(const Sparse<T1> &, const Sparse<T2> &);

        template<typename T1, typename T2>
        friend Sparse<decltype(T1() * T2())>
        operator-(const Sparse<T1> &, const Sparse<T2> &);

        template<typename T1, typename T2>
        friend Sparse<decltype(T1() * T2())>
        operator*(Sparse<T1> &, Sparse<T2> &);

        template<typename T1, typename T2>
        friend Sparse<decltype(T1() * T2())>
        operator/(const Sparse<T1> &, const Sparse<T2> &);

        template<typename T1, typename T2>
        friend Sparse<decltype(T1() * T2())>
        operator%(const Sparse<T1> &, const Sparse<T2> &);

        template<typename T1, typename T2>
        friend Sparse<decltype(T1() * T2())>
        operator+=(const Sparse<T1> &, const Sparse<T2> &);

        template<typename T1, typename T2>
        friend Sparse<decltype(T1() * T2())>
        operator-=(const Sparse<T1> &, const Sparse<T2> &);

        template<typename T1, typename T2>
        friend Sparse<decltype(T1() * T2())>
        operator*=(Sparse<T1> &, Sparse<T2> &);

        template<typename T1, typename T2>
        friend Sparse<decltype(T1() * T2())>
        operator/=(const Sparse<T1> &, const Sparse<T2> &);

        template<typename T1, typename T2>
        friend Sparse<decltype(T1() * T2())>
        operator%=(const Sparse<T1> &, const Sparse<T2> &);
    };


/*
 * ! Zeros are stored too, maybe we can strip them immediately here,
 * without additional complexity?
*/
    template<typename T>
    Sparse<T>::Sparse(const std::vector<T> &coefficients) {
        int index = coefficients.size();
        n = index - 1;
        std::ranges::for_each(coefficients | std::views::reverse, [&index, this](const T &p) {
            index--;
            if (p != 0) { monomials.push_back(Monomial<T>(p, index)); }
        });
    }

    template<typename T>
    Sparse<T>::Sparse(const std::vector<T> &&coefficients) {
        int index = coefficients.size();
        n = index - 1;
        std::ranges::for_each(coefficients | std::views::reverse, [&index, this](const T &p) {
            index--;
            if (p != 0) { monomials.push_back(Monomial<T>(std::move(p), index)); }
        });
    }

    template<typename T>
    Sparse<T>::Sparse(const std::list<T> &coefficients) {
        int index = coefficients.size();
        n = index - 1;
        std::ranges::for_each(coefficients | std::views::reverse, [&index, this](const T &p) {
            index--;
            if (p != 0) { monomials.push_back(Monomial<T>(p, index)); }
        });
        this->adjust();
    }

    template<typename T>
    Sparse<T>::Sparse(const std::list<T> &&coefficients) {
        int index = coefficients.size();
        n = index - 1;
        std::ranges::for_each(coefficients | std::views::reverse, [&index, this](const T &p) {
            index--;
            if (p != 0) { monomials.push_back(Monomial<T>(p, index)); }
        });
        this->adjust();
    }

    template<typename T>
    Sparse<T>::Sparse(const std::map<int, T, std::greater<>> &coefficients) {
        std::ranges::transform(coefficients, std::back_inserter(this->monomials),
                               [this](const std::pair<int, T> &p) {
                                   n = std::max(n, p.first);
                                   return Monomial<T>(p.second, p.first);
                               });
        this->adjust();
    }

    template<typename T>
    Sparse<T>::Sparse(const std::map<int, T, std::greater<>> &&coefficients) {
        std::ranges::transform(coefficients, std::back_inserter(this->monomials),
                               [this](std::pair<int, T> &p) {
                                   n = std::max(n, p.first);
                                   return Monomial<T>(std::move(p.second), p.first);
                               });
        this->adjust();
    }

    template<typename T>
    Sparse<T>::Sparse(const std::map<int, T> &coefficients) {
        std::ranges::transform(coefficients, std::front_inserter(this->monomials),
                               [this](const std::pair<int, T> &p) {
                                   n = std::max(n, p.first);
                                   return Monomial<T>(p.second, p.first);
                               });
        this->adjust();
    }

    template<typename T>
    Sparse<T>::Sparse(const std::map<int, T> &&coefficients) {
        std::ranges::transform(coefficients.begin(), coefficients.end(), std::front_inserter(this->monomials),
                               [this](std::pair<int, T> &p) {
                                   n = std::max(n, p.first);
                                   return Monomial<T>(std::move(p.second), p.first);
                               });
        this->adjust();
    }

    template<typename T>
    Sparse<T>::Sparse(const std::unordered_map<int, T> &&coefficients) {
        std::ranges::transform(coefficients, std::back_inserter(this->monomials),
                               [this](std::pair<int, T> &p) {
                                   n = std::max(n, p.first);
                                   return Monomial<T>(std::move(p.second), p.first);
                               });
        this->reorder();
        this->adjust();
    }

    template<typename T>
    Sparse<T>::Sparse(const std::unordered_map<int, T> &coefficients) {
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
        this->adjust();
    }


    template<typename T>
    Sparse<T>::Sparse(const std::vector<std::pair<int, T>> &coefficients) {
        /*
         * The input vector is assumed ordered by degree
         */
        std::ranges::transform(coefficients, std::back_inserter(this->monomials),
                               [this](const std::pair<int, T> &p) {
                                   n = std::max(n, p.first);
                                   return Monomial<T>(p.second, p.first);
                               });
        this->adjust();
    }

    template<typename T>
    Sparse<T>::Sparse(const std::vector<std::pair<int, T>> &&coefficients) {
        /*
         * The input vector is assumed ordered by degree
         */
        std::ranges::transform(coefficients, std::back_inserter(this->monomials),
                               [this](std::pair<int, T> &p) {
                                   n = std::max(n, p.first);
                                   return Monomial<T>(std::move(p.second), p.first);
                               });
        this->adjust();
    }

    template<typename T>
    Sparse<T>::Sparse(const std::list<std::pair<int, T>> &coefficients) {
        /*
         * The input list is assumed ordered by degree
         */
        std::ranges::transform(coefficients, std::back_inserter(this->monomials),
                               [this](const std::pair<int, T> &p) {
                                   n = std::max(n, p.first);
                                   return Monomial<T>(p.second, p.first);
                               });
        this->adjust();
    }

    template<typename T>
    Sparse<T>::Sparse(const std::list<std::pair<int, T>> &&coefficients) {
        /*
         * The input list is assumed ordered by degree
         */
        std::ranges::transform(coefficients, std::back_inserter(this->monomials), [this](std::pair<int, T> &p) {
            n = std::max(n, p.first);
            return Monomial<T>(std::move(p.second), p.first);
        });
        this->adjust();
    }

    template<typename T>
    Sparse<T>::Sparse(const Polynomial::Dense<T> &P, bool check_sparsity) {
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
    void Sparse<T>::derivative(int k) {
        if (k == 0) { return; }
        if (degree() < k) {
            *this = Sparse<T>();
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
    Sparse<T> derivative(Sparse<T> P, int k) {
        P.derivative(k);
        return P;
    }

    template<typename T>
    template<typename U>
    T Sparse<T>::operator()(const U &x) const {
        T result = 0;
        for (const auto &m: monomials) {
            result += m.coeff() * std::pow(x, m.degree());
        }
        return result;
    }

    template<typename T>
    std::ostream &operator<<(std::ostream &out, const Sparse<T> &P) {
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
    bool operator==(const Sparse<T1> &P, const Sparse<T2> &Q) {
        if (P.degree() != Q.degree()) { return false; }
        auto it1 = P.begin();
        auto it2 = Q.begin();
        while (it1 != P.end() && it2 != Q.end()) {
            if (*it1 != *it2) { return false; }
            ++it1;
            ++it2;
        }
//        return it1 == P.end() && it2 == Q.end();
        return true;
    }

    template<typename T1, typename T2>
    Sparse<decltype(T1() * T2())>
    operator+(const Sparse<T1> &lhs, const Sparse<T2> &rhs) {
        Sparse<decltype(T1() * T2())> result;
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
    Sparse<decltype(T1() * T2())>
    operator-(const Sparse<T1> &lhs, const Sparse<T2> &rhs) {
        Sparse<decltype(T1() * T2())> result(rhs);
        std::ranges::transform(result.monomials, begin(result.monomials),
                               [](const Monomial<decltype(T1() * T2())> &m) { return negate(m); });
        return lhs + result;
    }

    template<typename T1, typename T2>
    Sparse<decltype(T1() * T2())> operator*(Sparse<T1> &P, Sparse<T2> &Q) {
        std::priority_queue<Monomial<decltype(T1() * T2())>> heap;
        Sparse<decltype(T1() * T2())> result;
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
    std::pair<Sparse<T>, Sparse<T>>
    euclid_div(const Sparse<T> &P, const Sparse<T> &Q) {
        Sparse<T> quotient;
        auto remainder = P;
        while (!remainder.monomials.empty() && remainder.monomials.back().degree() >= Q.monomials.back().degree()) {

        }

        return std::pair<Sparse<T>, Sparse<T>>(quotient, remainder);
    }

    template<typename T1, typename T2>
    std::pair<Sparse<decltype(T1() * T2())>, Sparse<decltype(T1() * T2())>>
    heap_div(const Sparse<T1> &dividend, const Sparse<T2> &divisor) {
        Sparse<decltype(T1() * T2())> quotient;
        Sparse<decltype(T1() * T2())> remainder = dividend;
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
    Sparse<decltype(T1() * T2())> operator/(const Sparse<T1> &P, const Sparse<T2> &Q) {
        Sparse<decltype(T1() * T2())> quotient;
        auto remainder = P;

        while (!remainder.monomials.empty() && remainder.monomials.back().degree() >= Q.monomials.back().degree()) {
            auto leadCoeff = remainder.monomials.back().coeff() / Q.monomials.back().coeff();
            int leadDegree = remainder.monomials.back().degree() - Q.monomials.back().degree();
            quotient.monomials.push_back(Monomial<decltype(T1() * T2())>{leadCoeff, leadDegree});

            Sparse<decltype(T1() * T2())> product;
            for (const auto &mono: Q) {
                product.monomials.push_back(
                        Monomial<decltype(T1() * T2())>{mono.coeff() * leadCoeff, mono.degree() + leadDegree});
            }

            remainder = remainder - product;
        }
        return quotient;
    }

    template<typename T1, typename T2>
    Sparse<decltype(T1() * T2())> operator%(const Sparse<T1> &, const Sparse<T2> &) {
        return Sparse<decltype(T1() * T2())>();
    }

    template<typename T1, typename T2>
    Sparse<decltype(T1() * T2())> operator+=(const Sparse<T1> &A, const Sparse<T2> &B) {
        return A + B;
    }

    template<typename T1, typename T2>
    Sparse<decltype(T1() * T2())> operator-=(const Sparse<T1> &A, const Sparse<T2> &B) {
        return A - B;
    }

    template<typename T1, typename T2>
    Sparse<decltype(T1() * T2())> operator*=(Sparse<T1> &A, Sparse<T2> &B) {
        return A * B;
    }

    template<typename T1, typename T2>
    Sparse<decltype(T1() * T2())>
    operator/=(const Sparse<T1> &, const Sparse<T2> &) {/* TODO */}

    template<typename T1, typename T2>
    Sparse<decltype(T1() * T2())>
    operator%=(const Sparse<T1> &, const Sparse<T2> &) {/* TODO */}

}

#endif //POLYNOMIALSHPC_SPARSEPOLYNOMIALS_HPP
