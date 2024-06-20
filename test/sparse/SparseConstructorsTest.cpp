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
        auto zero2 = Monomial(0, 0);
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
        std::map<int, double, std::greater<>> mg;
        mg[4] = 5;
        mg[3] = 1;
        mg[2] = 0;
        mg[1] = 3;
        mg[0] = 6.3;
        Sparse MG(mg);
        std::unordered_map<int, double> um;
        um[4] = 5;
        um[3] = 1;
        um[2] = 0;
        um[1] = 3;
        um[0] = 6.3;
        Sparse UM(um);
        std::vector<std::pair<int, double>> pv({{4, 5},
                                                {3, 1},
                                                {2, 0},
                                                {1, 3},
                                                {0, 6.3}});
        Sparse PV(pv);
        std::list<std::pair<int, double>> pl({{4, 5},
                                              {3, 1},
                                              {2, 0},
                                              {1, 3},
                                              {0, 6.3}});
        Sparse PL(pl);
        EXPECT_EQ(V, L);
        EXPECT_EQ(L, M);
        EXPECT_EQ(M, MG);
        EXPECT_EQ(MG, UM);
        EXPECT_EQ(UM, PV);
        EXPECT_EQ(PV, PL);

    }
}

