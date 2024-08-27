//
// Created by Emmanuel Leclercq on 13/04/2024.
//
#include <iostream>
#include "DensePolynomials.hpp"
#include "gtest/gtest.h"

using namespace Polynomial;
namespace {
    TEST(SolverTest, DegreeOne) {
        std::vector<int> v1 = {1, 6};
        Dense p1(v1, true);
        std::cout << p1 << std::endl;
        std::cout << "p[0]=" << p1[0] << " p[1]=" << p1[1] << " p[2]=" << p1[2] << std::endl;
        Dense<double> p2({2, 3}, true);

        std::vector<int> roots_p1 = solveRoots<int, int>(p1);
        EXPECT_EQ(roots_p1, v1);
        std::vector<double> roots_p2 = solveRoots<double, double>(p2);


    }
}
