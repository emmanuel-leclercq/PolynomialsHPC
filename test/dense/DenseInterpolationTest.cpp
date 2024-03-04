//
// Created by Emmanuel Leclercq on 27/02/2024.
//
#include <iostream>
#include "Polynomials.hpp"


#include "gtest/gtest.h"

namespace {
    TEST(InterpolationTest, Basic) {
        auto interpolation = interpolate<double>({{0.0, 0.0},
                                                  {2.0, 4.0},
                                                  {4.0, 16.0},
                                                  });
        Polynomial<double> P({0, 0, 1});
        EXPECT_EQ(interpolation, P);
    }
}