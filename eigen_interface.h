//
// CopyRight Chong Peng 2017
//

#ifndef CODE_EXAMPLE_EIGEN_INTERFACE_H_
#define CODE_EXAMPLE_EIGEN_INTERFACE_H_

#include "type_def.h"

namespace code_example {

/**
 * Interface of Eigen Matrix to use DavidsonDiag classes
 *
 * \@tparam T the element type in Eigen Matrix
 *
 */

template <typename T>
RowMatrix<T> copy_zero(const RowMatrix<T>& a){
    return RowMatrix<T>::Zero(a.rows(), a.cols());
}

template <typename T>
void scale(RowMatrix<T>& y, T a) {
  y = a * y;
}

template <typename T>
void axpy(RowMatrix<T>& y, T a, const RowMatrix<T>& x){
  y = y + a*x;
}

template <typename T>
T dot_product(const RowMatrix<T>& x, const RowMatrix<T>& y){
  T dot = x.cwiseProduct(y).sum();
  return dot;
}

template <typename T>
T norm2(const RowMatrix<T>& a){
  return std::sqrt(a.squaredNorm());
}

} // namespace code_example

#endif // CODE_EXAMPLE_EIGEN_INTERFACE_H_
