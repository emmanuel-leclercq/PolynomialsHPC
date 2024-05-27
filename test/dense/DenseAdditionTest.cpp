//
// Created by Emmanuel Leclercq on 26/12/2023.
//
#include "Polynomials.hpp"

#include "gtest/gtest.h"

using namespace Polynomial;
namespace {
    TEST(AdditionTest, Trivial) {
        Dense<int> p1{{6, 3, 0, 1, 5}};
        Dense<int> p2{{1, 0, 1}};

        Dense<int> p3{{7, 3, 1, 1, 5}};

        EXPECT_EQ(p1 + p2, p3);
    }
}
