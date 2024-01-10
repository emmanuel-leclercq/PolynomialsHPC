//
// Created by Emmanuel Leclercq on 26/12/2023.
//
#include <iostream>
#include "Polynomials.hpp"

#include <gtest/gtest.h>
namespace {
    TEST(DerivativeTest, Trivial) {
    Polynomial p1({6, 3, 0, 1, 5});
    Polynomial p2({3, 0, 3, 20});
    Polynomial p3({0,3,0,1,5});
    Polynomial p4({0,0,0, 0.25,1});

    EXPECT_EQ(p1.derivative(), p2);
    EXPECT_EQ(p1.derivative(-1),p3);
    EXPECT_EQ(derivative(p1),p2);
    EXPECT_EQ(derivative(p1,-1),p4);
    EXPECT_EQ(derivative(p1,6),Polynomial<double>());

    }
}