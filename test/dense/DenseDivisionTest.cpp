//
// Created by Emmanuel Leclercq on 26/12/2023.
//
#include <gtest/gtest.h>
#include <iostream>
#include "Polynomials.hpp"

#include "gtest/gtest.h"
namespace {
    TEST(DivisionTest, Quotient) {
    Polynomial<int> p1{{7, 3, 1, 1, 5}};
    Polynomial<int> p2{{1, 0, 1}};
    Polynomial<int> quotient{{5, 5, 5}};

    EXPECT_EQ(p1 / p2, quotient);
    }

    TEST(DivisionTest,Remainder){
        Polynomial<int> p1{{7, 3, 1, 1, 5}};
        Polynomial<int> p2{{1, 0, 1}};
        Polynomial<int> remainder{{2, 5}};

        EXPECT_EQ(p1 % p2, remainder);
    }
}