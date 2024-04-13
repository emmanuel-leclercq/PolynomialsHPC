//
// Created by Emmanuel Leclercq on 13/04/2024.
//
#include <iostream>
#include "Polynomials.hpp"
#include "gtest/gtest.h"

namespace {
    TEST(SolverTest, DegreeOne) {
        Polynomial<int> p1{{6, 1}};
        Polynomial<double> p2({2, 3}, true);

        std::vector<double> Rootsp1 = solveRoots<int, double>(p1);
        std::vector<double> Rootsp2 = solveRoots<double, double>(p2);


    }
}
