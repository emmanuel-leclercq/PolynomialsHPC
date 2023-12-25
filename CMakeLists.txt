cmake_minimum_required(VERSION 3.26)
project(PolynomialsHPC)

set(CMAKE_CXX_STANDARD 20)
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address")

# add_library(${PROJECT_NAME} SHARED Polynomials.cpp)

#add_executable(Dense TestPolynomial.cpp)

add_executable(Sparse TestSparsePolynomial.cpp)

#add_executable(bench multiplicationBenchmark.cpp)

