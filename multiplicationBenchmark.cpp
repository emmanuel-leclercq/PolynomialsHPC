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
    std::cout << p1 * p2 << std::endl;
    std::cout << (ans == (p1 * p2)) << std::endl;
    for (int i = 0; i < 50; i++) {
        auto P = generateRandomIntPolynomial(5, 0, 1);
        auto Q = generateRandomIntPolynomial(3, 0, 1);
        if ((karatsuba(P, Q)) != P * Q) {
            std::cout << "FFT: " << fftmultiply(P, Q) << std::endl;
            std::cout << "brute force=" << P * Q << std::endl;
            std::cout << "P=" << P << std::endl;
            std::cout << "Q=" << Q << std::endl;
            std::cout << "P.degree()=" << P.degree() << std::endl;
            std::cout << "Q.degree()=" << Q.degree() << std::endl << std::endl;

        }
    }

    Timer timer;
    auto R = generateRandomIntPolynomial(4096, -100, 100);
    auto S = generateRandomIntPolynomial(4096, -100, 100);
    timer.start();
    karatsuba(R, S);
    timer.finish();
    std::cout << "Polynomial deg 1000 multiplication using karatsuba: " << timer() << "ms" << std::endl;
    timer.start();
    R * S;
    timer.finish();
    std::cout << "Polynomial deg 1000 multiplication using brute force: " << timer() << "ms" << std::endl;




//    auto P = generateRandomIntPolynomial(1000, -100, 100);
//    auto Q = generateRandomIntPolynomial(1000, -100, 100);

    return 0;
}