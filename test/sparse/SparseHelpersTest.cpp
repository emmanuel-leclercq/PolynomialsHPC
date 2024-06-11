//
// Created by Emmanuel_Leclercq on 11/06/2024.
//
#include "SparsePolynomials.hpp"

#include "gtest/gtest.h"
#include <list>

using namespace Polynomial;
namespace {
    TEST(HelpersTest, Trivial) {
        //reorder
        Sparse<int> p1(std::list<int>{6, 3, 0, 1, 5});
        /*
        Sparse<int> p2(std::vector<std::pair<int, int>>{
                {0, 3},
                {1, 1},
                {3, 1},
                {5, 1},
                {6, 1}});

        p1.reorder();
        p2.reorder();
        EXPECT_EQ(p1, p2);
         */
    }
}