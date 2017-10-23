//
// Created by Chong Peng on 10/20/17.
//

#ifndef CODE_EXAMPLE_STD_VECTOR_INTERFACE_H_
#define CODE_EXAMPLE_STD_VECTOR_INTERFACE_H_

#include <vector>

typedef int32_t integer4;

/// BLAS function declarations

extern "C" {
void dscal_(const integer4*, const double*, double*, const integer4*);
void daxpy_(const integer4 *, const double *, const double *, const integer4 *,
            double *, const integer4 *);
double ddot_(const integer4 *, const double *, const integer4 *, const double *,
             const integer4 *);
double dnrm2_(const integer4* , const double *, const integer4*);
void dgemv_(const char*, const integer4*, const integer4*, const double*,
            const double*, const integer4*, const double*, const integer4*,
            const double*, double*, const integer4*);
};

namespace code_example {

/**
 * Interface of std::vector to use DavidsonDiag classes
 *
 * currently only support double as element type
 */

inline std::vector<double> copy_and_zero(const std::vector<double> &a) {
  auto result = std::vector<double>(a.size(), 0.0);
  result.shrink_to_fit();
  return result;
}

inline void scale(std::vector<double> &y, double a) {
  int32_t y_size = y.size();
  int32_t inc = 1;
  dscal_(&y_size, &a, y.data(), &inc);
}

inline void axpy(std::vector<double> &y, double a,
                 const std::vector<double> &x) {
  int32_t y_size = y.size();
  assert(y_size == x.size());
  int32_t inc = 1;
  daxpy_(&y_size, &a, x.data(), &inc, y.data(), &inc);
}

inline double dot_product(const std::vector<double> &x,
                          const std::vector<double> &y) {
  int32_t x_size = x.size();
  assert(x_size == y.size());
  int32_t inc = 1;
  return ddot_(&x_size, x.data(), &inc, y.data(), &inc);
}

inline double norm2(const std::vector<double> &a) {
  int32_t a_size = a.size();
  int32_t inc = 1;
  return dnrm2_(&a_size, a.data(), &inc);
}

} // namespace code_example

#endif // CODE_EXAMPLE_STD_VECTOR_INTERFACE_H_
