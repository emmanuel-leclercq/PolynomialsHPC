//
// Created by Emmanuel Leclercq on 26/12/2023.
//
#include <iostream>
#include "SparsePolynomials.hpp"

#include "gtest/gtest.h"

namespace {

    TEST(MonomialTest, basic) {
        Polynomial::Monomial<int> m;
        Polynomial::Monomial m1(2, 2);
        Polynomial::Monomial m2(1, 3);
        auto zero = Polynomial::Monomial<int>(0);


        EXPECT_EQ(m, 0);
        EXPECT_EQ(zero, 0);
        EXPECT_TRUE(m1 < m2);
    }
}

