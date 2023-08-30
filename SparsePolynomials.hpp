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

template <typename T>
class Monomial;

template <typename T>
class SparsePolynomial;

template <typename T>
bool operator<(const Monomial<T> &u, const Monomial<T> &v) { return u.n < v.n; };

template <typename T>
bool operator==(const Monomial<T> &u, const Monomial<T> &v) { return u.n == v.n && u.coefficient == v.coefficient; };

template <typename T>
bool operator==(const Monomial<T> &u, const T &t) { return u.coefficient == t; };

template <typename T>
std::ostream &operator<<(std::ostream &, const Monomial<T> &);

template <typename T>
std::ostream &operator<<(std::ostream &, const SparsePolynomial<T> &);

template <typename T>
class Monomial
{
private:
    int n;
    T coefficient;

public:
    Monomial() : coefficient(0), n(-1) {}

    Monomial(const T &a, const int &m = 0) : coefficient(a), n(m) {}

    Monomial(const T &&a, const int &m = 0) : coefficient(a), n(m) {}

    /* The goal is to order by degree
     * ! Do not auto reformat this with VSCode or erro
     */
    friend bool operator< <>(const Monomial<T> &u, const Monomial<T> &v);

    friend bool operator==<>(const Monomial<T> &u, const Monomial<T> &v);

    friend bool operator==<>(const Monomial<T> &u, const T &t);

    friend std::ostream &operator<<<>(std::ostream &, const Monomial<T> &);
};

template <typename T>
class SparsePolynomial
{
private:
    std::list<Monomial<T>> monomials;
    int n;
    bool is_sorted=false;

public:
    /*
     * We provide a handful of constructors for convenience and flexibility,
     * considering certain data types make more sense for a sparse implementation,
     * and to reduce the need for preprocessing existing variables
     */
    SparsePolynomial() : monomials(Monomial<T>()), n(-1), is_sorted(true) {}

    /*
     * Constructors assuming input is not 'sparse', which spears us
     * the need to know the degree beforehand
     */
    explicit SparsePolynomial(const std::list<T> &coefficients);

    explicit SparsePolynomial(const std::list<T> &&coefficients);

    explicit SparsePolynomial(const std::vector<T> &coefficients);

    explicit SparsePolynomial(const std::vector<T> &&coefficients);

    /* It makes more sense to use degrees as keys considering the
     * vectorial structural of polynomials spaces
     */
    explicit SparsePolynomial(const std::vector<std::pair<int, T>> &coefficients);

    explicit SparsePolynomial(const std::vector<std::pair<int, T>> &&coefficients);

    explicit SparsePolynomial(const std::list<std::pair<int, T>> &coefficients);

    explicit SparsePolynomial(const std::list<std::pair<int, T>> &&coefficients);

    explicit SparsePolynomial(const std::map<int, T> &coefficients);

    explicit SparsePolynomial(const std::map<int, T> &&coefficients);

    explicit SparsePolynomial(const std::unordered_map<int, T> &coefficients);

    explicit SparsePolynomial(const std::unordered_map<int, T> &&coefficients);

    /*
     * Conversion from dense to sparse, with optional check for sparsity (more than half zeros)
     */
    explicit SparsePolynomial(const Polynomial<T> &P, bool check_sparsity = false);

    //  Default destructor will be used, standard data types
    ~SparsePolynomial() = default;

    [[nodiscard]] int degree() const { return n; }

    [[nodiscard]] Monomial<T> dominant() const { return *monomials.begin(); }

    [[nodiscard]] Monomial<T> &dominant() { return *monomials.begin(); }

    void reorder()
    {
        monomials.sort([](Monomial<T> &P, Monomial<T> &Q)
                       { return Q < P; });
        is_sorted = true;
    }

    void adjust()
    {
        monomials.remove_if([](const Monomial<T> &M)
                            { return M == 0; });
    };

    friend std::ostream &operator<<<>(std::ostream &, const SparsePolynomial<T> &);
};

template <typename T>
SparsePolynomial<T>::SparsePolynomial(const std::map<int, T> &coefficients)
{
    std::transform(coefficients.begin(), coefficients.end(), std::back_inserter(this->monomials),
                   [](const std::pair<int, T> &p)
                   { return Monomial<T>(p.second, p.first); });
}

template <typename T>
SparsePolynomial<T>::SparsePolynomial(const std::map<int, T> &&coefficients)
{
    std::transform(coefficients.begin(), coefficients.end(), std::back_inserter(this->monomials),
                   [](std::pair<int, T> &p)
                   { return Monomial<T>(std::move(p.second), p.first); });
    delete coefficients;
}

template <typename T>
SparsePolynomial<T>::SparsePolynomial(const std::unordered_map<int, T> &&coefficients)
{
    std::transform(coefficients.begin(), coefficients.end(), std::back_inserter(this->monomials),
                   [](std::pair<int, T> &p)
                   { return Monomial<T>(std::move(p.second), p.first); });
    delete coefficients;
}

template <typename T>
SparsePolynomial<T>::SparsePolynomial(const std::unordered_map<int, T> &coefficients)
{
    std::transform(coefficients.begin(), coefficients.end(), std::back_inserter(this->monomials),
                   [](const std::pair<int, T> &p)
                   { return Monomial<T>(p.second, p.first); });
}

/*
 * !Problem with the zeros
*/
template <typename T>
SparsePolynomial<T>::SparsePolynomial(const std::vector<T> &coefficients)
{
    int index = coefficients.size();
    n = index;
    std::transform(coefficients.begin(), coefficients.end(), std::back_inserter(this->monomials),
                   [&index](const T &p)
                   { return Monomial<T>(p, index--); });
}

template <typename T>
SparsePolynomial<T>::SparsePolynomial(const std::vector<T> &&coefficients)
{
    int index = coefficients.size();
    n = index;
    std::transform(coefficients.begin(), coefficients.end(), std::back_inserter(this->monomials),
                   [&index](T &p)
                   { return Monomial<T>(std::move(p), index--); });
}

template <typename T>
SparsePolynomial<T>::SparsePolynomial(const std::list<T> &coefficients) {}

template <typename T>
SparsePolynomial<T>::SparsePolynomial(const std::list<T> &&coefficients) {}

template <typename T>
SparsePolynomial<T>::SparsePolynomial(const std::vector<std::pair<int, T>> &coefficients) {}

template <typename T>
SparsePolynomial<T>::SparsePolynomial(const Polynomial<T> &P, bool check_sparsity)
{
    if (check_sparsity)
    {
        int zeros = 0;
        for (auto it = P.coefficients.begin(); it != P.coefficients.end(); ++it)
        {
            if (is_zero(*it))
            {
                zeros++;
            }
        }
        if (zeros > P.coefficients.size() / 2)
        {
            throw std::invalid_argument("The polynomial is not sparse enough");
        }
    }
    int index = P.coefficients.size();
    n = index;
    for (auto &c : P.coefficients)
    {
        if (!is_zero(c))
        {
            monomials.emplace_back(c, index);
        }
        index--;
    }
}

template <typename T>
std::ostream &operator<<(std::ostream &out, const Monomial<T> &p)
{

    if (p.n == 0)
    {
        out << p.coefficient;
    }
    else if (is_zero(p.coefficient))
    {
        out << "0";
    }
    else if (is_one(p.coefficient))
    {
        out << "x^" << p.n;
    }
    else
    {
        out << p.coefficient << "x^" << p.n;
    }
    return out;
}

template <typename T>
std::ostream &operator<<(std::ostream &out, const SparsePolynomial<T> &P)
{
    for (auto it = P.monomials.begin(); it != P.monomials.end(); ++it)
    {
        out << *it;
        if (it != P.monomials.end() && std::next(it) != P.monomials.end())
        {
            out << " + ";
        }
    }
    {
    }
    return out;
}

#endif //POLYNOMIALSHPC_SPARSEPOLYNOMIALS_HPP
