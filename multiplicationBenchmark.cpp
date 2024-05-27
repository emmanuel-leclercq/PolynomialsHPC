/*
 * The goal here is to use either Karatsuba, FFT or brute force
 * multiplication whenever they are faster, depending on the input degrees
 */
#include "Polynomials.hpp"
#include "utils.hpp"
#include <iostream>
#include <vector>

using namespace Polynomial;
int main() {

    std::vector<int> v1{6, 3, 0, 1, 5};
    Dense<int> p1(v1);

    std::vector<int> v2{1, 0, 1};
    Dense<int> p2(v2);
    auto ans = fftmultiply(p1, p2);
    auto cant = cantor(p1, p2);
    std::cout << p1 << std::endl;
    std::cout << p2 << std::endl;
    std::cout << ans << std::endl;
    std::cout << cant << std::endl;
    std::cout << p1 * p2 << std::endl;
    std::cout << (ans == (p1 * p2)) << std::endl;
    for (int i = 0; i < 50; ++i) {
        auto R = generateRandomIntPolynomial(31, -1000, 1000); //above deg 31 cantor fails??
        auto S = generateRandomIntPolynomial(100, -1000, 1000);
        std::cout << (R * S == cantor(R, S)) << std::endl;

    }

    Timer timer;


    auto R = generateRandomIntPolynomial(10000, -1000, 1000);
    auto S = generateRandomIntPolynomial(10000, -1000, 1000);

    timer.start();
    R * S;
    timer.finish();
    std::cout << "Polynomial deg 10000 multiplication using brute force: " << timer() << "ms" << std::endl;
    timer.start();
    fftmultiply(R, S);
    timer.finish();
    std::cout << "Polynomial deg 10000 multiplication using FFT: " << timer() << "ms" << std::endl;
//    std::cout << fftmultiply(R, S);


//    auto P = generateRandomIntPolynomial(1000, -100, 100);
//    auto Q = generateRandomIntPolynomial(1000, -100, 100);

    return 0;
}