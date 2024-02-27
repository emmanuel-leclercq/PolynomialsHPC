//
// Created by Emmanuel Leclercq on 26/12/2023.
//
#include <gtest/gtest.h>
#include <iostream>
#include "Polynomials.hpp"

#include "gtest/gtest.h"
namespace {
    TEST(MultiplicationTest, Zero) {
        Polynomial<int> p1{{6, 3, 0, 1, 5}};
        Polynomial<int> p2;
        EXPECT_EQ(p1*p2,p2);
    }
    TEST(MultiplicationTest, Normal) {
        Polynomial<int> p1{{6, 3, 0, 1, 5}};
        Polynomial<int> p2{{1, 0, 1}};
        Polynomial<int> p3{{6, 3, 6, 4, 5, 1, 5}};
        EXPECT_EQ(p1*p2,p3);
    }
}