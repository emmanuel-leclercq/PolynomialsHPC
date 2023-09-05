//
// Created by Emmanuel Leclercq on 05/09/2023.
//
/*
 * We try a different sparse implementation using a map.
 * The point being the ability to have monomials natively ordered by degree.
 * We will try to replicqte the performances of the linked list implementation
 * both in terms of memory and speed.
 */
#ifndef POLYNOMIALSHPC_MAPSPARSE_HPP
#define POLYNOMIALSHPC_MAPSPARSE_HPP

#include <map>

template<typename T>
class SparseTest {
private:
    std::map<int, T> monomials;
public:
    explicit SparseTest(std::map<int,T> m): monomials(m){};
};

#endif //POLYNOMIALSHPC_MAPSPARSE_HPP
