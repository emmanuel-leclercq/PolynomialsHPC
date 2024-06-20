//
// Created by Emmanuel Leclercq on 26/12/2023.
//
#include "DensePolynomials.hpp"

#include "gtest/gtest.h"

using namespace Polynomial;
namespace {
    TEST(MultiplicationTest, Zero) {
        Dense<int> p1{{6, 3, 0, 1, 5}};
        Dense<int> p2;
        EXPECT_EQ(p1 * p2, p2);
    }

    TEST(MultiplicationTest, Normal) {
        Dense<int> p1{{6, 3, 0, 1, 5}};
        Dense<int> p2{{1, 0, 1}};
        Dense<int> p3{{6, 3, 6, 4, 5, 1, 5}};
        EXPECT_EQ(p1 * p2, p3);
    }
}