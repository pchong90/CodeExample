# CodeExample
CodeExample is a repository for my example codes

#### Generic Solver of Davidson Algorithm 
The solver works for symmetric and non-symmetric system with real eigenvalues.  
The generic design makes it possible to use this solver with different data structure.   
Two examples are provided that use `Eigen::Matrix`(eigen_test.cpp) and `std::vector`(std_vector_test.cpp) as the type of the eigenvectors.

Dependency: Eigen, BLAS, CMake  
To compile: `cmake .`, then `make eigen_test` or `make std_vector_test`
