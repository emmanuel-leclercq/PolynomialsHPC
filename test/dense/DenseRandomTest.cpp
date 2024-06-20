//
// Created by Emmanuel Leclercq on 26/12/2023.
//
#include "DensePolynomials.hpp"

#include "gtest/gtest.h"

using namespace Polynomial;
namespace {
    TEST(RandomTest, Unif) {
        auto P = generateRandomIntPolynomial(100, 0, 1);
        /*
         * This test is not ideal: there is a non-zero
         * (although very low) probability that the polynomial is zero
         */
        EXPECT_NE(P.degree(), 0);
    }

    TEST(RandomTest, CustomGenerator) {
        std::random_device rd;
        std::default_random_engine def_generator(rd());
        std::normal_distribution double_distribution(0.0, 1.0);

        auto Q = generateRandomPolynomial<double>(30, double_distribution, def_generator);

        EXPECT_EQ(Q.degree(), 30);
    }
}