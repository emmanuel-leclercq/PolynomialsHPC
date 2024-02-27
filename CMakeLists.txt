cmake_minimum_required(VERSION 3.20)
project(PolynomialsHPC)

set(CMAKE_CXX_STANDARD 20 REQUIRED)
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -v")

include(FetchContent)
FetchContent_Declare(
        googletest
        URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
)
FetchContent_MakeAvailable(googletest)

include_directories(/src)
include_directories(/test)
include_directories(${CMAKE_SOURCE_DIR})

# add_library(${PROJECT_NAME} SHARED Polynomials.cpp)

#add_executable(Dense TestPolynomial.cpp)

#add_executable(Sparse TestSparsePolynomial.cpp)

#add_executable(bench multiplicationBenchmark.cpp)

enable_testing()
file(GLOB_RECURSE SOURCES "test/dense/*.cpp" )
add_executable(tests ${SOURCES})
#
target_link_libraries(tests gtest_main)
include(GoogleTest)
gtest_discover_tests(tests)


