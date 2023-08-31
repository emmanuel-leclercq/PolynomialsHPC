#include "Polynomials.hpp"
#include <iostream>
#include <vector>
#include <chrono>

using std::cout;
using std::endl;
using std::vector;

int main() {

    /*
    * Testing constructors
    */

    Polynomial<double> q;
    cout << "Default polynomial degree : " << q.degree() << endl;

    vector<int> v1{6, 3, 0, 1, 5};
    Polynomial<int> p1(v1);

    vector<int> v2{1, 0, 1};
    Polynomial<int> p2(v2);

    std::complex<double> a(2.0, 1.0);
    std::complex<double> b(0.0, 1.0);
    vector<std::complex<double>> vc{a, b};

    Polynomial<std::complex<double>> pc(vc);

    std::complex<int> one = 1;
    cout << "Is one one : " << is_one(one) << endl;

    cout << is_zero(a) << endl;
    cout << Polynomial(a, 3) << endl;
    cout << endl;

    cout << "p1 : " << p1 << endl;
    cout << "p2 : " << p2 << endl;

    /*
    * Testing basic operations
    */
    Polynomial<int> sum = p1 + p2;
    Polynomial<int> diff = p1 - p2;
    Polynomial<int> diff2 = p2 - p1;
    Polynomial<int> prod = p1 * p2;
    cout << "p1+p2 : " << sum << endl;
    cout << "p1-p2 : " << diff << endl;
    cout << "p2-p1 : " << diff2 << endl;
    cout << "p1*p2 : " << prod << endl;

    // Division and remainder
    Polynomial<int> div = p1 / p2;
    Polynomial<int> reste = p1 % p2;
    cout << "Quotient : " << div << endl;
    cout << "Reste : " << reste << endl;

    // functiòn evaluation
    cout << "p1(2) : " << p1(2) << endl
         << "p2(3) " << p2(3) << endl;
    // Random polynomials generation
    auto P = generateRandomPolynomial(100, 0, 1);
    auto Q = generateRandomPolynomial(100, -50, 50);

    cout << "Is random polynomial P sparse? " << P.is_sparse() << endl;
    cout << "P: " << P << endl;
    auto start = std::chrono::high_resolution_clock::now();
    cout << "Dominant coefficient of random polynomial product: " << (P * Q).dominant() << endl;

    auto end = std::chrono::high_resolution_clock::now();
    cout << "time it takes for the product (and cout): ";
    cout << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() << "ms" << endl;

    return 0;
}
