#include "SparsePolynomials.hpp"
#include "Polynomials.hpp"
#include <iostream>
#include <vector>
#include <map>
#include <unordered_map>

using std::cout;
using std::endl;
using std::list;
using std::vector;

int main() {
    /*
     * testing Monomial implementation
    */
    Monomial<int> m1(1, 2);
    Monomial<int> m2(1, 3);
    cout << "m1: " << m1 << endl;
    cout << "m1<m2 in terms of degree? " << (m1 < m2) << endl;
    auto zero = Monomial<int>(0);
    cout << "m1=0? " << (m1 == 0) << endl;
    cout << "monomial zero=0? " << (zero == 0) << endl;

    /*
    * testing Sparse Polynomial implementation
    */
    cout << "constructing V from vector" << endl;
    vector<double> vector1{6.3, 3, 0, 1, 5};
    SparsePolynomial<double> V(vector1, false);
    cout << "V: " << V << endl;
    V.adjust();
    cout << "V adjusted: " << V << endl;
    cout << "V.degree(): " << V.degree() << endl;
    cout << "Dominant monomial: " << V.dominant() << endl;

    cout << "constructing L from list<int>" << endl;
    list<int> list1{6, 3, 0, 1, 5};
    SparsePolynomial<int> L(list1, false);
    cout << "L: " << L << endl;

    cout << "constructing M from map<int,double>" << endl;
    std::map<int, double> map1;
    map1[17] = 15.3;
    map1[175] = -5;
    map1[0] = 777;
    SparsePolynomial<double> M(map1, false);
    cout << "M: " << M << endl;

    cout << "constructing U from unordered_map<int,double>" << endl;
    std::unordered_map<int, double> umap1;
    umap1[17] = 1555.3;
    umap1[175] = -5555;
    umap1[0] = 777;
    SparsePolynomial<double> U(umap1, false);
    cout << "U: " << U << endl;

    cout << "constructing U from Polynomial<double>" << endl;
    Polynomial<double> poly1(vector1);
    SparsePolynomial<double> P(poly1, false);
    cout << "U: " << P << endl;

    return 0;
}