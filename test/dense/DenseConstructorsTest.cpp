//
// Created by Emmanuel Leclercq on 26/12/2023.
//
#include <iostream>
#include "Polynomials.hpp"

#include "gtest/gtest.h"
namespace {
    TEST(ConstructorTest, Trivial) {
    Polynomial<int> p1{{0, 0, 0, 0, 5}};
    Polynomial<int> p2(5, 4);
    Polynomial<double> p3{{1.5, 2.5}, true};
    Polynomial<double> p4{{3.75, -4.0, 1}};

    EXPECT_EQ(p1, p2);
    EXPECT_EQ(p3, p4);
}
}