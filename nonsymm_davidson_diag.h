//
// Copyright Chong Peng 2017
//

#ifndef CODE_EXAMPLE_NONSYMM_DAVIDSON_DIAG_H_
#define CODE_EXAMPLE_NONSYMM_DAVIDSON_DIAG_H_

#include "symm_davidson_diag.h"

namespace code_example {

namespace detail {

template <typename T, typename Vector> struct EigenPair {
  using element_type = T;
  using result_type = Vector;

  element_type eigen_value;
  result_type eigen_vector;

  /// constructor
  EigenPair(const element_type &value, const result_type &vector)
      : eigen_value(value), eigen_vector(vector) {}

  /// move constructor
  EigenPair(const element_type &&value, const result_type &&vector)
      : eigen_value(std::move(value)), eigen_vector(std::move(vector)) {}

  ~EigenPair() = default;

  // sort by eigen value
  bool operator<(const EigenPair &other) const {
    return eigen_value < other.eigen_value;
  }
};

} // namespace detail

// clang-format off
/**
*  Davidson Algorithm for Non-symmetric systems with real eigenvalues
*
*  this Class inherits from SymmDavidsonDiag class and overrides two
*  differences between NonSymmDavidsonDiag and SymmDavidsonDiag
*
*  - void compute_new_subspace(): how the subspace is computed
*  - void eigen_solve_subspace(result_type &E, RowMatrix<element_type> &C): how the subspace is diagnolized
*
*/
//clang-format on
template <typename Array>
class NonSymmDavidsonDiag : public SymmDavidsonDiag<Array> {

public:
  using typename SymmDavidsonDiag<Array>::element_type;
  using typename SymmDavidsonDiag<Array>::result_type;
  using typename SymmDavidsonDiag<Array>::value_type;

  /**
   *  Constructor See SymmDavidsonDiag for documentation
   */
  NonSymmDavidsonDiag(unsigned int n_roots, unsigned int n_guess = 2,
                      unsigned int max_n_guess = 4)
      : SymmDavidsonDiag<Array>(n_roots, n_guess, max_n_guess) {}

  using SymmDavidsonDiag<Array>::solve;
  using SymmDavidsonDiag<Array>::extrapolate;

private:

  /// this override function to compute_new_subspace for non-symmetrix systems
  void compute_new_subspace() override {

    // size of original subspace
    const auto n_s = this->subspace_.cols();
    // size of new subspace
    const auto n_v = this->B_.size();
    // size of new vector
    const auto n_b = n_v - n_s;

    // compute new subspace
    // G will be replicated Eigen Matrix
    {
      RowMatrix<element_type> G = RowMatrix<element_type>::Zero(n_v, n_v);
      // reuse stored subspace
      G.block(0, 0, n_s, n_s) << this->subspace_;
      // initialize new value
      for (std::size_t i = 0; i < n_b; ++i) {
        const auto ii = i + n_s;
        for (std::size_t j = 0; j <= ii; ++j) {
          G(ii, j) = dot_product(this->B_[ii], this->HB_[j]);
          if (ii != j) {
            G(j, ii) = dot_product(this->B_[j], this->HB_[ii]);
          }
        }
      }

      this->subspace_ = G;
    }
  }

  /// this override function to do eigen_solve_subspace for non-symmetrix systems
  void eigen_solve_subspace(result_type &E,
                            RowMatrix<element_type> &C) override {
    // nonsymmetric matrix
    Eigen::EigenSolver<RowMatrix<element_type>> es(this->subspace_);

    // sort eigen values
    std::vector<detail::EigenPair<element_type, result_type> > eg;
    {
      RowMatrix<element_type> v = es.eigenvectors().real();
      EigenVector<element_type> e = es.eigenvalues().real();

      if (es.info() != Eigen::Success) {
        throw std::runtime_error("Eigen::EigenSolver Failed!\n");
      }

      const std::size_t n_v = v.cols();
      for (std::size_t i = 0; i < n_v; ++i) {
        eg.emplace_back(e[i], v.col(i));
      }

      std::sort(eg.begin(), eg.end());
    }

    // obtain final eigen value and eigen vector
    for (std::size_t i = 0; i < this->n_roots_; ++i) {
      E[i] = eg[i].eigen_value;
      C.col(i) = eg[i].eigen_vector;
    }
  }
};

} // namespace code_example

#endif // CODE_EXAMPLE_NONSYMM_DAVIDSON_DIAG_H_
