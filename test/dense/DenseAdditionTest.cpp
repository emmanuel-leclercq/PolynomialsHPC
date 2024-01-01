//
// Created by Emmanuel Leclercq on 26/12/2023.
//
#include <iostream>
#include "Polynomials.hpp"


#include "gtest/gtest.h"
namespace {
    TEST(AdditionTest, Trivial) {
        Polynomial<int> p1{{6, 3, 0, 1, 5}};
        Polynomial<int> p2{{1, 0, 1}};

        Polynomial<int> p3{{7, 3, 1, 1, 5}};

        EXPECT_EQ(p1 + p2, p3);
    }
}
