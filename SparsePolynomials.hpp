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

template <typename T>
class monomial;

template <typename T>
bool operator<(const monomial<T> &u, const monomial<T> &v) { return u.n < v.n; };

template <typename T>
class monomial
{
private:
    int n;
    T coefficient;

public:
    monomial() : coefficient(0), n(-1) {}

    monomial(const T &a, const int &m = 0) : coefficient(a), n(m) {}

    monomial(const T &&a, const int &m = 0) : coefficient(a), n(m) {}

    //  Ordering by degree
    friend bool operator<<T>(const monomial<T> &u, const monomial<T> &v);
};

template <typename T>
class SparsePolynomial
{
private:
    std::list<monomial<T>> monomials;
    int n{};

public:
    /*
     * We provide a handful of constructors for convenience and flexibility,
     * considering certain data types make more sense for a sparse implementation,
     * and to reduce the need for preprocessing existing variables
     */
    SparsePolynomial() : monomials(monomial<T>()), n(-1) {}

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

    explicit SparsePolynomial(const std::list<std::pair<int, T>> &coefficients);

    explicit SparsePolynomial(const std::vector<std::pair<int, T>> &&coefficients);

    explicit SparsePolynomial(const std::list<std::pair<int, T>> &&coefficients);

    explicit SparsePolynomial(const std::map<int, T> &coefficients);

    explicit SparsePolynomial(const std::map<int, T> &&coefficients);

    explicit SparsePolynomial(const std::unordered_map<int, T> &coefficients);

    explicit SparsePolynomial(const std::unordered_map<int, T> &&coefficients);

    //  Default constructors will be used, standard data types
    ~SparsePolynomial() = default;
};

template <typename T>
SparsePolynomial<T>::SparsePolynomial(const std::map<int, T> &coefficients)
{
    std::transform(coefficients.begin(), coefficients.end(), std::back_inserter(this->monomials),
                   [](const std::pair<int, T> &p)
                   { return monomial(p.second, p.first); });
}

template <typename T>
SparsePolynomial<T>::SparsePolynomial(const std::map<int, T> &&coefficients)
{
    std::transform(coefficients.begin(), coefficients.end(), std::back_inserter(this->monomials),
                   [](std::pair<int, T> &p)
                   { return monomial(std::move(p.second), p.first); });
    delete coefficients;
}

template <typename T>
SparsePolynomial<T>::SparsePolynomial(const std::unordered_map<int, T> &&coefficients)
{
    std::transform(coefficients.begin(), coefficients.end(), std::back_inserter(this->monomials),
                   [](std::pair<int, T> &p)
                   { return monomial(std::move(p.second), p.first); });
    delete coefficients;
}

template <typename T>
SparsePolynomial<T>::SparsePolynomial(const std::unordered_map<int, T> &coefficients)
{
    std::transform(coefficients.begin(), coefficients.end(), std::back_inserter(this->monomials),
                   [](const std::pair<int, T> &p)
                   { return monomial(p.second, p.first); });
}

template <typename T>
SparsePolynomial<T>::SparsePolynomial(const std::vector<T> &coefficients)
{
    int index = coefficients.size();
    std::transform(coefficients.begin(), coefficients.end(), std::back_inserter(this->monomials),
                   [](const T &p)
                   { return monomial(p, index); });
}

template <typename T>
SparsePolynomial<T>::SparsePolynomial(const std::vector<T> &&coefficients)
{
    int index = coefficients.size();
    std::transform(coefficients.begin(), coefficients.end(), std::back_inserter(this->monomials),
                   [](T &p)
                   { return monomial(std::move(p), index); });
}

template <typename T>
SparsePolynomial<T>::SparsePolynomial(const std::list<T> &coefficients) {}

template <typename T>
SparsePolynomial<T>::SparsePolynomial(const std::list<T> &&coefficients) {}

template <typename T>
SparsePolynomial<T>::SparsePolynomial(const std::vector<std::pair<int, T>> &coefficients) {}

#endif //POLYNOMIALSHPC_SPARSEPOLYNOMIALS_HPP
