//
// Created by Emmanuel Leclercq on 26/12/2023.
//
#include "Polynomials.hpp"

#include "gtest/gtest.h"

using namespace Polynomial;
namespace {
    TEST(ConstructorTest, Trivial) {
    Dense<int> p1{{0, 0, 0, 0, 5}};
    Dense<int> p2(5, 4);
    Dense<double> p3{{1.5, 2.5}, true};
    Dense<double> p4{{3.75, -4.0, 1}};

    EXPECT_EQ(p1, p2);
    EXPECT_EQ(p3, p4);
}
}