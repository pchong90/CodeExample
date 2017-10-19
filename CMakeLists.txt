cmake_minimum_required(VERSION 3.0.0)

project(CodeExample)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


# find Eigen3
FIND_PACKAGE(Eigen3 REQUIRED)
INCLUDE_DIRECTORIES(${EIGEN3_INCLUDE_DIR})

# find TBB
#FIND_PACKAGE(TBB REQUIRED)


set(src davidson_diag.h gram_schmidt.h)


add_executable(eigen_test EXCLUDE_FROM_ALL eigen_test.cpp ${src})