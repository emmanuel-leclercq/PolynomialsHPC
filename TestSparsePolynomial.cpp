// TestSparsePolynomial.cpp

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
    Monomial<int> m;
    Monomial<int> m1(2, 2);
    Monomial<int> m2(1, 3);
    cout << "m1: " << m1 << endl;
    cout << "m1<m2 in terms of degree? " << (m1 < m2) << endl;
    auto zero = Monomial<int>(0);
    cout << "m1=0? " << (m1 == 0) << endl;
    cout << "monomial zero=0? " << (zero == 0) << endl;
    cout << "default monomial: " << m << endl;

    /*
    * testing Sparse Polynomial implementation
    */
    SparsePolynomial<int> basic;
    cout << "default sparse polynomial" << basic << endl;

    cout << "constructing V from vector" << endl;
    vector<double> vector1{6.3, 3, 0, 1, 5};
    SparsePolynomial<double> V(vector1);
    cout << "V: " << V << endl;
    V.adjust();
    cout << "V adjusted: " << V << endl;
    cout << "V.degree(): " << V.degree() << endl;
    cout << "Dominant monomial: " << V.dominant() << endl;

    cout << "constructing L from list<int>" << endl;
    list<int> list1{6, -3, 0, 1, 5};
    SparsePolynomial<int> L(list1);
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
    cout << "dense version: " << poly1 << endl;
    SparsePolynomial<double> P(poly1, false);
    cout << "Sparse version P: " << P << endl;
    cout << "Testing copy 2nd derivative: P'=" << derivative<double>(P, 2) << endl;
    P.derivative();
    cout << "Testing derivative P=" << P << endl;
    P.derivative(-1);
    cout << "Testing antiderivative P=" << P << endl;
    P.derivative(3);
    cout << "Testing 3rd derivative P=" << P << endl;
    P.derivative(-4);
    cout << "Testing 4th antiderivative P=" << P << endl;

    SparsePolynomial<int> R;
    R.add({1, 5});
    R.add({1, 2});
    R.add({3, 0});
    cout << "R: " << R << endl;
    SparsePolynomial<int> S;
    S.add({3, 3});
    S.add({-1, 0});
    cout << "S: " << S << endl;
    cout << "R*S: " << R * S << endl;

    cout << "Testing addition M+P=" << M + P << endl;
    cout << "Testing addition P+M=" << P + M << endl;
    cout << "Testing subtraction M-P=" << M - P << endl;
    cout << "Testing subtraction P-M=" << P - M << endl;
    cout << "Testing multiplication P*M=" << P * M << endl;
//    cout << "Testing division M/P=" << M / P << endl;
//    cout << "Testing division M%P=" << M % P << endl;

    //profiling addition
//    auto A_dense=generateRandomIntPolynomial(100, 0, 1);
//    auto A_sparse=SparsePolynomial<int>(A_dense, false);
//    auto B_dense=generateRandomIntPolynomial(100, 0, 1);
//    auto B_sparse=SparsePolynomial<int>(B_dense, false);
//    Timer timer;
//    timer.start();
//    A_sparse+B_sparse;
//    timer.finish();
//    cout<<"classic addition time: "<<timer()<<endl;
//    timer.start();
//    heap_plus<int>(A_sparse, B_sparse);
//    timer.finish();
//    cout<<"heap addition time: "<<timer()<<endl;
//
//    unsigned long long brute=0, heap=0;
//    for (int i = 0; i < 100; ++i) {
//        //average time of addition brute force versus heap:
//        auto A_dense=generateRandomIntPolynomial(1000, 0, 1);
//        auto A_sparse=SparsePolynomial<int>(A_dense, false);
//        auto B_dense=generateRandomIntPolynomial(1000, 0, 1);
//        auto B_sparse=SparsePolynomial<int>(B_dense, false);
//        Timer timer;
//        timer.start();
//        A_sparse+B_sparse;
//        timer.finish();
//        brute+=timer();
//        timer.start();
//        heap_plus<int>(A_sparse, B_sparse);
//        timer.finish();
//        heap+=timer();
//    }
//    cout<<"average time of addition brute force versus heap: "<<brute/100<<" "<<heap/100<<endl;

    auto Q = generateRandomIntPolynomial(100, 0, 1);
    cout << "Q and its degree: " << Q << ", " << Q.degree() << endl;
    cout << "Is Q sparse? " << Q.is_sparse() << endl;
    SparsePolynomial<int> Q_sparse(Q, true);
    cout << "Q as sparse and its degree: " << Q_sparse << ", " << Q_sparse.degree() << endl;

    return 0;
}