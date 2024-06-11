//
// Created by Emmanuel Leclercq on 26/12/2023.
//
#include "SparsePolynomials.hpp"

#include "gtest/gtest.h"
#include <list>

using namespace Polynomial;
namespace {
    TEST(DerivativeTest, Trivial) {
        Sparse<int> p1(std::list<int>{6, 3, 0, 1, 5});
        Sparse<double> p1Double(std::list<double>{6, 3, 0, 1, 5});
        Sparse<int> p2(std::list<int>{3, 0, 3, 20});
        Sparse<int> p3(std::list<int>{0, 3, 0, 1, 5});
        Sparse<double> p4(std::list<double>{0, 6, 1.5, 0, 0.25, 1});

        p1.derivative();
        EXPECT_EQ(p1, p2);
        p1.derivative(-1);
        EXPECT_EQ(p1, p3);

        EXPECT_EQ(derivative(p1), p2);
        EXPECT_EQ(derivative(p1Double, -1), p4);
        EXPECT_EQ(derivative(p1, 6), Sparse<int>());

    }
}