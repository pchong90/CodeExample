# CodeExample
CodeExample is a repository for my example codes

#### Generic Solver of Davidson Algorithm 
Dependency: Eigen3, BLAS, CMake  
To compile: `cmake .`, then `make eigen_test` or `make std_vector_test`


This Davidson Algorithm solver works for symmetric(`symm_davidson_diag.h`) and non-symmetric(`nonsymm_davidson_diag.h`) system with real eigenvalues. The generic design makes it possible to use this solver with eigenvectors from different data structure classes. The user needs to provide the following interface functions, where the Array is the template type which stands for the type of eigenvector  
 - `element_type_traits<Array>::type` defines the `element_type` in Array (for example float or double)
 - `Array copy_and_zero(const Array& a)`  make a new Array based on a, while initialize all elements as zero
 - `element_type dot_product(const Array& a, const Array& b)` return dot product of Array a and b
 - `void scale(Array& y , element_type a)` scale Array y by factor a
 - `void axpy(Array&y , element_tye a, const Array& x)` performs y += a*x
 - `element_type norm2(const Array& a)` return the L^2 norm of Array a


Two examples are provided that use `Eigen::Matrix`(`eigen_test.cpp`) and `std::vector`(`std_vector_test.cpp`) as the type of the eigenvectors. The `element_type_traits` functions are defined in `type_traits.h`. All other interface functions are defined in `eigen_interface.h` and `std_vector_interface.h`.


