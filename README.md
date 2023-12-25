# PolynomialsHPC

This is an attempt to build a headers-only high performance single-variable polynomial functions libraries,
over any ring (the ring doesn't have to be compatible with the Fourier transform). The dense polynomial
class serves as a benchmark for the sparse polynomial class, which aims to achieve better memory performances.

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
   the latter. Two implementations -> P.derivative() and derivative(P) for convenience (the latter doesn't modify P)
9. Root finder (TODO)

### Sparse Polynomial class

A polynomial is considered sparse when more than half its coefficients are zeros.
Unlike the dense polynomial class, we store coefficients in a linked list (x^999+1 is stored as two pairs
(coefficient, degree) in a list instead of a 1000-sized vector).
Assuming t is the number of elements in the list representing a sparse polynomial, our goal is to approach O(t)
complexity for most operations, which is not always possible. In fact, some operations are simply not possible,
such as factorization (except for specific cases).

#### Operations available

1. Addition
2. Subtraction
3. Multiplication (Johnson, Monagan, Pearce heap algorithm)
4. Division (TODO)
5. Evaluation (TODO)
6. nth derivatives
7. Root finder (TODO)

#### Resources

Modern Computer Algebra, Gathen & Gerhard