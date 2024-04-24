//
// Created by Emmanuel Leclercq on 13/04/2024.
//
#include <iostream>
#include "Polynomials.hpp"
#include "gtest/gtest.h"

namespace {
    TEST(SolverTest, DegreeOne) {
        std::vector<int> v1 = {1, 6};
        Polynomial<int> p1(v1, true);
        std::cout << p1 << std::endl;
        std::cout<<"p[0]="<<p1[0]<<" p[1]="<<p1[1]<<" p[2]="<<p1[2]<<std::endl;
        Polynomial<double> p2({2, 3}, true);

        std::vector<int> Rootsp1 = solveRoots<int, int>(p1);
        EXPECT_EQ(Rootsp1, v1);
//        std::vector<double> Rootsp2 = solveRoots<double, double>(p2);


    }
}
