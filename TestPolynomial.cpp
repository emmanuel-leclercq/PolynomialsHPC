// TestPolynomial.cpp

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

    vector<double> v3{6, 3, 0, 1, 5};
    Polynomial<double> p3(v3);

    std::complex<double> a(2.0, 1.0);
    std::complex<double> b(0.0, 1.0);
    vector<std::complex<double>> vc{a, b};
    Polynomial<std::complex<double>> pc(vc);

    std::complex<int> one = 1;
    cout << "Is one one : " << is_one(one) << endl;

    cout << "is 2+i zero : " << is_zero(a) << endl;
    cout << Polynomial(a, 3) << endl;
    cout << endl;

    cout << "p1 : " << p1 << endl;
    cout << "p2 : " << p2 << endl;
    cout << endl;

    Polynomial<double> withRoots({1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12}, true);
    cout << "Building polynomial from roots 1,2,3,...,12: " << withRoots << endl << endl;
    std::vector<double> points{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13};
    auto ans{withRoots.multipointEval(points)};
    cout << "checking correctness with multipoint evaluation: ";
    for (auto x: ans) { cout << x << " "; }
    cout << endl << endl;

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
    cout << endl;

    // Division and remainder
    Polynomial<int> div = p1 / p2;
    Polynomial<int> remainder = p1 % p2;
    cout << "Quotient : " << div << endl;
    cout << "Remainder : " << remainder << endl;
    cout << endl;
    // function evaluation
    cout << "p1(2) : " << p1(2) << endl;
    cout << "p2(3) : " << p2(3) << endl;
    cout << endl;

    //derivative and antiderivative
    p3.derivative();
    cout << "Test p1.derivative() = " << p3 << endl;
    cout << "p1.degree() = " << p3.degree() << endl;
    p3.derivative(-1);
    cout << "Test p1.derivative(-1) = " << p3 << endl;
    cout << "p1.degree() = " << p3.degree() << endl;
    cout << "test derivative(p1) = " << derivative<double>(p3) << endl;
    cout << "test derivative(p1,-1) = " << derivative<double>(p3, -1) << endl;
    cout << endl;
    cout << "Test interpolation (0,0), (1,1) (3,9): ";
    vector<std::pair<double, double>> data = {{0, 0},
                                              {2, 4},
                                              {4, 16}};
    auto interpolation = interpolate<double>({{0, 0},
                                              {2, 4},
                                              {4, 16}});
    cout << interpolation << endl;
    // Random polynomials generation
    auto P = generateRandomIntPolynomial(100, 0, 1);

    std::random_device rd;
    std::mt19937 mt_generator(rd());
    std::default_random_engine def_generator(rd());
    std::normal_distribution<double> double_distribution(0.0, 1.0);

    auto Q = generateRandomPolynomial<double>(30, double_distribution, def_generator);

    cout << "P : " << P << endl;
    cout << "Q : " << Q << endl;
    Timer timer;
    timer.start();
    auto dom = (P * Q).dominant();
    timer.finish();
    cout << "Dominant coefficient of random polynomial product : " << dom << endl;
    cout << "Time it takes for the product : ";
    cout << timer() << "ms" << endl;
    cout << endl;

    return 0;
}
