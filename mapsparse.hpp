//
// Created by Emmanuel Leclercq on 05/09/2023.
//
/*
 * We try a different sparse implementation using a map.
 * The point being the ability to have monomials natively ordered by degree.
 * We will try to replicate the performances of the linked list implementation
 * both in terms of memory and speed.
 */
#ifndef POLYNOMIALSHPC_MAPSPARSE_HPP
#define POLYNOMIALSHPC_MAPSPARSE_HPP

#include <map>

template<typename T>
class SparseTest {
public:
    std::map<int, T, std::greater<int>> monomials;

    explicit SparseTest(std::map<int, T> m) {
        for (auto [key, val]: m) { monomials[key] = val; }
    };

    explicit SparseTest(std::map<int, T, std::greater<int>> m) : monomials(m) {};

    int degree() { return *this->begin()->first; }
};

#endif //POLYNOMIALSHPC_MAPSPARSE_HPP
