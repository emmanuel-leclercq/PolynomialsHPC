//
// Created by Emmanuel Leclercq on 27/02/2024.
//
#include "Polynomials.hpp"

#include "gtest/gtest.h"

using namespace Polynomial;
namespace {
    TEST(InterpolationTest, Basic) {
        auto interpolation = interpolate<double>({{0.0, 0.0},
                                                  {2.0, 4.0},
                                                  {4.0, 16.0},
                                                  });
        Dense<double> P({0, 0, 1});
        EXPECT_EQ(interpolation, P);
    }
}