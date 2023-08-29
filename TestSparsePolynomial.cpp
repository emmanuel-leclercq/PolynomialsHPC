#include "SparsePolynomials.hpp"
#include <iostream>
#include <vector>

using std::cout;
using std::endl;
using std::vector;

int main()
{
    vector<int> v1{6, 3, 0, 1, 5};
    SparsePolynomial<int> p1(v1);
    return 0;
}