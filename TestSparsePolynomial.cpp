#include "SparsePolynomials.hpp"
#include <iostream>
#include <vector>

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
    cout << "constructing p1 from vector" << endl;
    vector<int> v1{6, 3, 0, 1, 5};
    SparsePolynomial<int> p1(v1, false);
    cout << "p1: " << p1 << endl;
    p1.adjust();
    cout << "p1 adjusted: " << p1 << endl;
    cout << "p1.degree(): " << p1.degree() << endl;
    cout << "Dominant monomial: " << p1.dominant() << endl;

    return 0;
}