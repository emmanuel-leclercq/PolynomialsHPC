//
// Created by Emmanuel Leclercq on 26/12/2023.
//
#include "DensePolynomials.hpp"

#include "gtest/gtest.h"

using namespace Polynomial;
namespace {
    TEST(DivisionTest, Quotient) {
        Dense<int> p1{{7, 3, 1, 1, 5}};
        Dense<int> p2{{1, 0, 1}};
        Dense<int> quotient{{5, 5, 5}};

        EXPECT_EQ(p1 / p2, quotient);
    }

    TEST(DivisionTest, Remainder) {
        Dense<int> p1{{7, 3, 1, 1, 5}};
        Dense<int> p2{{1, 0, 1}};
        Dense<int> remainder{{2, 5}};

        EXPECT_EQ(p1 % p2, remainder);
    }
}