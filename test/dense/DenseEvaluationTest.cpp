//
// Created by Emmanuel Leclercq on 26/12/2023.
//
#include <iostream>
#include "Polynomials.hpp"

#include "gtest/gtest.h"
namespace {
    TEST(EvaluationTest, SinglePoint) {
        Polynomial<int> P({1, 2, 3, 4,}, true);
        Polynomial<double> Q(1.5, 2);
        EXPECT_EQ(P(1), 0);
        EXPECT_EQ(Q(2.0), 6.0);
    }


    TEST(EvaluationTest, MultiPoint) {
        std::vector<int> roots{1, 2, 3, 4,};
        Polynomial<int> P(roots, true);
        std::vector<int> ans{0,0,0,0};
        EXPECT_EQ(P(roots),ans);
    }
}