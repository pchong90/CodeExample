//
// Created by Chong Peng on 10/19/17.
//

#ifndef CODE_EXAMPLE_TYPE_TRAITS_H_
#define CODE_EXAMPLE_TYPE_TRAITS_H_

#include <vector>
#include "type_def.h"

namespace code_example {

template <typename T, typename Enabler = void>
struct element_type_traits {};

template <typename T>
struct element_type_traits<RowMatrix<T>> {
  typedef T type;
};

template <typename T>
struct element_type_traits<std::vector<T>> {
  typedef T type;
};

} // namespace code_example

#endif // CODE_EXAMPLE_TYPE_TRAITS_H_
