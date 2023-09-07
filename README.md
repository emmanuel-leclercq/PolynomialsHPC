# PolynomialsHPC

This is an attempt to build a headers-only high performance polynomial functions libraries, single parameter for now.

### Polynomial (dense) class

This class represents polynomial coefficients as a vector, degrees are implied.
We make heavy use of the <algorithm> and <numeric> libraries for efficiency.
The adjust() method removes trailing zeros, basic accessors are provided.
the is_sorted parameter is left to user to set for efficiency.
Addition between non-ordered sparse polynomials results in a sorted sparse polynomial (no sorting overhead).

#### Operations available :

1. addition: O(min(n,m)) implementation, we take the largest (in terms of degree) of lhs and rhs, then add the terms
2. Subtraction: similar method to addition, but not using addition for speed
3. Multiplication: implementation in O(n*m) for now, will try FFT for O(nlog(n))
4. Division: euclid algorithm
5. Evaluation: Horner's method
6. Factorization (TODO)
7. Interpolation: takes a vector<std::pair<T x,T y>> as input and returns a dense polynomial<T>
8. nth derivatives: includes antiderivative (n<0), both derivative and antiderivative are O(n), we use a swap trick for
   the latter. Two implementations -> P.derivative() and derivative(P) for convenience (the latter doesnt modify P)
9. Root finder (TODO)
10. Plotting (TODO)

### Sparse Polynomial class

A polynomial is considered sparse when more than half its coefficients are zeros.
Unlike the dense polynomial class, we store coefficients in a linked list (x^999+1 is stored in a 2-elements lists
instead of a 1000 elements vector).
Assuming t is the number of elements in the list representing a sparse polynomial, our goal is to approach O(t)
complexity for most operations, which is not always possible. We detail below the implementations (TODO)

#### Operations available

1. Addition (TODO)
2. Subtraction (TODO)
3. Multiplication (TODO)
4. Division (TODO)
5. Evaluation (TODO)
6. Factorization (TODO)
7. Interpolation (TODO)
8. nth derivatives (TODO)
9. Root finder (TODO)
10. Plotting (TODO)
