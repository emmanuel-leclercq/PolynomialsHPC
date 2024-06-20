//
// Created by Emmanuel Leclercq on 26/12/2023.
//
#include <iostream>

#include <vector>
#include <list>
#include <map>
#include <unordered_map>

#include "SparsePolynomials.hpp"
#include "DensePolynomials.hpp"

#include "gtest/gtest.h"

using namespace Polynomial;
namespace {

    TEST(ConstructorTest, monomials) {
        Monomial<int> M;
        auto zero1 = Monomial(0);
        auto zero2 = Monomial(0.0, 0);
        EXPECT_EQ(M, zero1);
        EXPECT_EQ(zero1, zero2);
    }

    TEST(ConstructorTest, Sparse) {
        Sparse<int> basic;
        std::vector<double> v{6.3, 3, 0, 1, 5};
        Sparse V(v);
        std::list<double> l{6.3, 3, 0, 1, 5};
        Sparse L(l);
        std::map<int, double> m;
        m[4] = 5;
        m[3] = 1;
        m[2] = 0;
        m[1] = 3;
        m[0] = 6.3;
        Sparse M(m);

        EXPECT_EQ(V, L);
        EXPECT_EQ(L, M);

    }
}

