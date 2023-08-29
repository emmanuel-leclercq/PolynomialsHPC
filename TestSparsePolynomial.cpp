#include "SparsePolynomials.hpp"
#include <iostream>
#include <vector>

using std::cout;
using std::endl;
using std::list;
using std::vector;

int main()
{
    /*
     * testing Monomial implementation
    */
    Monomial<int> m1(1, 2);
    Monomial<int> m2(1, 3);
    cout << "m1: " << m1 << endl;
    cout << "m1<m2 ? " << (m1 < m2) << endl;
    auto zero = Monomial<int>(0);
    cout << "m1=0? " << (m1 == 0) << endl;
    cout << "zero=0? " << (zero == 0) << endl;

    vector<int> v1{6, 3, 0, 1, 5};
    SparsePolynomial<int> p1(v1);
    cout << "p1: " << p1 << endl;
    p1.reorder();
    cout << "p1 reordered: " << p1 << endl;
    p1.adjust();
    cout << "p1 adjusted: " << p1 << endl;
    cout << "p1.degree(): " << p1.degree() << endl;
    cout << "Dominant monomial: " << p1.dominant() << endl;

    return 0;
}