//
// CopyRight Chong Peng 2017
//

#ifndef CODE_EXAMPLE_TYPE_DEF_H_
#define CODE_EXAMPLE_TYPE_DEF_H_

#include <Eigen/Dense>

/// typedef for EigenVector
template <typename T>
using EigenVector = ::Eigen::Matrix<T, 1, ::Eigen::Dynamic>;

/// typedef for RowMatrix
template <typename T>
using RowMatrix =
::Eigen::Matrix<T, ::Eigen::Dynamic, ::Eigen::Dynamic, ::Eigen::RowMajor>;

#endif //CODE_EXAMPLE_TYPE_DEF_H_
