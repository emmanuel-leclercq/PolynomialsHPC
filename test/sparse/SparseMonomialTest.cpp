//
// Created by Emmanuel Leclercq on 26/12/2023.
//
#include <iostream>
#include "SparsePolynomials.hpp"

#include "gtest/gtest.h"

using namespace Polynomial;
namespace {

    TEST(MonomialTest, basic) {
        Monomial<int> m;
        Monomial m1(2, 2);
        Monomial m2(1, 3);
        auto zero = Monomial<int>(0);


        EXPECT_EQ(m, 0);
        EXPECT_EQ(zero, 0);
        EXPECT_TRUE(m1 < m2);
    }
}

