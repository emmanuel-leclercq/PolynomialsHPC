cmake_minimum_required(VERSION 3.26)
project(PolynomialsHPC)

set(CMAKE_CXX_STANDARD 20)

# add_library(${PROJECT_NAME} SHARED Polynomials.cpp)

add_executable(${PROJECT_NAME} *.cpp)