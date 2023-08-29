#include "Polynomials.hpp"
#include <iostream>
#include <vector>
#include <chrono>

using std::cout;
using std::endl;
using std::vector;

int main()
{
    auto start = std::chrono::high_resolution_clock::now();
    //Première partie
    Polynomial<double> q;
    cout << "Degre du polynome construit par defaut : " << q.degree() << endl;

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

    //Deuxième partie
    // Somme, différénce, produit
    Polynomial<int> sum = p1 + p2;
    Polynomial<int> diff = p1 - p2;
    Polynomial<int> diff2 = p2 - p1;
    Polynomial<int> prod = p1 * p2;
    cout << "Somme : " << sum << endl;
    cout << "Différence p1-p2 : " << diff << endl;
    cout << "Différence p2-p1 : " << diff2 << endl;
    cout << "Produit : " << prod << endl;

    // Division et reste
    Polynomial<int> div = p1 / p2;
    Polynomial<int> reste = p1 % p2;
    cout << "Quotient : " << div << endl;
    cout << "Reste : " << reste << endl;

    // Evaluation en un point
    cout << "p1(2) : " << p1(2) << endl
         << "p2(3) " << p2(3) << endl;

    auto end = std::chrono::high_resolution_clock::now();
    cout << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() << "ms" << endl;

    return 0;
}