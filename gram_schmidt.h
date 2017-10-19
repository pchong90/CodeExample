//
// Copyright Chong Peng 2017
//

#ifndef CODE_EXAMLE_GRAM_SCHMIDT_H_
#define CODE_EXAMLE_GRAM_SCHMIDT_H_

#include <vector>

namespace code_example {

/**
 * Gram-Schmidt Algorithm
 *
 */


/**
 *   orthonomalize all Array in vector V, starting from position 0
 */
template <typename Array>
void gram_schmidt(std::vector<Array> &V, std::size_t start = 0) {
  const auto original_k = V.size();
  const auto k = V.size();
  std::size_t n_neglect = 0;

  assert(start < k);

  for (std::size_t i = start; i < k; ++i) {
    for (std::size_t j = 0; j < i; ++j) {
      const auto tmp = dot_product(V[i], V[j]);
      axpy(V[i], -tmp, V[j]);
    }

    // normalize
    auto norm = norm2(V[i]);
    scale(V[i], 1.0 / norm2(V[i]));
  }
}

/**
 *  vector V1 is already orthonormalized
 *  orthonormalize V2 with respect to V1
 */
template <typename Array>
void gram_schmidt(const std::vector<Array> &V1, std::vector<Array> &V2) {
  const auto original_k = V2.size();
  const auto k1 = V1.size();
  auto k2 = V2.size();
  std::size_t n_neglect = 0;

  for (std::size_t i = 0; i < k2; ++i) {
    // loop over all vector in V1
    for (std::size_t j = 0; j < k1; ++j) {
      auto tmp = dot_product(V2[i], V1[j]);
      axpy(V2[i], -tmp, V1[j]);
    }

    // loop over other vector in V2
    for (std::size_t k = 0; k < i; ++k) {
      auto tmp = dot_product(V2[i], V2[k]);
      axpy(V2[i], -tmp, V2[k]);
    }

    auto norm = norm2(V2[i]);
    scale(V2[i], 1.0 / norm);
  }
}

} // namespace code_example
#endif // CODE_EXAMLE_GRAM_SCHMIDT_H_
