/*
 * The goal here is to use either Karatsuba, FFT or brute force
 * multiplication whenever they are faster, depending on the input degrees
 */
#include "Polynomials.hpp"
#include "utils.hpp"
#include <iostream>
#include <vector>

int main() {

    std::vector<int> v1{6, 3, 0, 1, 5};
    Polynomial<int> p1(v1);

    std::vector<int> v2{1, 0, 1};
    Polynomial<int> p2(v2);
    auto ans = fftmultiply(p1, p2);
    std::cout << p1 << std::endl;
    std::cout << p2 << std::endl;
    std::cout << ans << std::endl;
    std::cout << p1 * p2;

//    auto P = generateRandomIntPolynomial(1000, -100, 100);
//    auto Q = generateRandomIntPolynomial(1000, -100, 100);

    return 0;
}