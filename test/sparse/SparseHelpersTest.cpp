//
// Created by Emmanuel_Leclercq on 11/06/2024.
//
#include "SparsePolynomials.hpp"

#include "gtest/gtest.h"
#include <list>
#include <vector>

using namespace Polynomial;
namespace {
    TEST(HelpersTest, Trivial) {
        //reorder
        Sparse p1(std::list<int>{6, 3, 0, 1, 5});

        std::vector<std::pair<int, int>> pv({{4, 5},
                                             {3, 1},
                                             {2, 0},
                                             {1, 3},
                                             {0, 6},});
        Sparse p2(pv);
        EXPECT_EQ(p1, p2); // implicit call to .reorder() and .adjust() removing null coefficients
        Monomial<int> b(8, 8);
        p1.add(b);
        p1.reorder();
        EXPECT_FALSE(p1 == p2);
        EXPECT_EQ(p1[8], b);


    }
}