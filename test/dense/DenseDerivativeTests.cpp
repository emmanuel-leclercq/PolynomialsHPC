//
// Created by Emmanuel Leclercq on 26/12/2023.
//
#include <iostream>
#include "Polynomials.hpp"

#include <gtest/gtest.h>

namespace {
    TEST(DerivativeTest, Trivial) {
        Polynomial<int> p1({6, 3, 0, 1, 5});
        Polynomial<double> p1Double({6, 3, 0, 1, 5});
        Polynomial<int> p2({3, 0, 3, 20});
        Polynomial<int> p3({0, 3, 0, 1, 5});
        Polynomial<double> p4({0, 6, 1.5, 0, 0.25, 1});

        p1.derivative();
        EXPECT_EQ(p1, p2);
        p1.derivative(-1);
        EXPECT_EQ(p1, p3);

        EXPECT_EQ(derivative(p1), p2);
        EXPECT_EQ(derivative(p1Double, -1), p4);
        EXPECT_EQ(derivative(p1, 6), Polynomial<int>());

    }
}