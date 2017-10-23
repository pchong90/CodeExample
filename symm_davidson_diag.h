//
// Copyright Chong Peng 2017
//

#ifndef CODE_EXAMPLE_DAVIDSON_DIAG_H_
#define CODE_EXAMPLE_DAVIDSON_DIAG_H_

#include <deque>
#include <stdexcept>
#include <vector>

#include "gram_schmidt.h"
#include "type_def.h"

/**
 *
 */

namespace code_example {

// clang-format off
/**
 * \brief An Example Code of Generic Davidson Algorithm Implementation For Symmetric System
 *
 * Solves eigen value problem <tt> Hx = ex </tt> for the lowest n eigen value and eigen vector
 *
 * it starts with orthogonal guess vector B {b1, b2, ... bn}
 * the eigen vector is a linear combination of B
 * extrapolate() will update the vector B and store new x
 *
 * \tparam Array array type, which is the type of the eigen vector and guess vector
 *
 * \c Array must provide the following interface functions:
 *
 * \c element_type_traits<Array>::type must be defined, which gives the element type in Array
 *
 *
 * - `Array copy_and_zero(const Array& a)`  make a new Array based on a, while initialize all elements as zero
 * - `element_type dot_product(const Array& a, const Array& b)` return dot product of 2 Array a and b
 * - `void scale(Array& y , element_type a)` scale Array y by factor a
 * - `void axpy(Array&y , element_tye a, const Array& x)` performs y += a*x
 * - `element_type norm2(const Array& a)` return the L^2 norm of Array a
 *
 */
// clang-format on


template <typename Array>
class SymmDavidsonDiag {
public:
  using element_type = typename element_type_traits<Array>::type;
  using result_type = EigenVector<element_type>;
  using value_type = std::vector<Array>;

public:
  /**
   *
   * @param n_roots number of lowest roots to solve
   *
   * @param n_guess number of eigen vector per root at subspace collapse,
   * default is 2
   *
   * @param max_n_guess max number of guess vector per root, default is 4
   *
   */
  SymmDavidsonDiag(unsigned int n_roots, unsigned int n_guess = 2,
                   unsigned int max_n_guess = 4)
      : n_roots_(n_roots), n_guess_(n_guess), max_n_guess_(max_n_guess),
        eigen_vector_(), HB_(), B_(), subspace_() {}

  virtual ~SymmDavidsonDiag() {
    eigen_vector_.clear();
    HB_.clear();
    B_.clear();
    subspace_.resize(0, 0);
  }

  /**
   *
   * @tparam Operator  operator that computes the product of H*B
   * @tparam Pred operator that precondition the residual
   *
   * @param guess initial guess vector
   * @param op    op(B) should compute HB
   * @param pred  preconditioner, pred(residucal) will precondition the residual
   * @param convergence   convergence threshold
   * @param max_iter  max number of iteration allowd
   *
   * @return solved eigen values
   */
  template <typename Operator, typename Pred>
  EigenVector<element_type> solve(value_type &guess, const Operator &op,
                                  const Pred &pred, double convergence,
                                  std::size_t max_iter) {
    double norm_e = 1.0;
    double norm_r = 1.0;
    std::size_t iter = 0;

    EigenVector<element_type> eig = EigenVector<element_type>::Zero(n_roots_);

    while (iter < max_iter && (norm_r > convergence || norm_e > convergence)) {

      // compute product of H with guess vector
      value_type HC = op(guess);

      EigenVector<element_type> eig_new, norms;
      std::tie(eig_new, norms) = extrapolate(HC, guess, pred);

      EigenVector<element_type> delta_e = (eig - eig_new);
      delta_e = delta_e.cwiseAbs();
      norm_e =
          *std::max_element(delta_e.data(), delta_e.data() + delta_e.size());
      norm_r = *std::max_element(norms.data(), norms.data() + norms.size());

//      std::cout << "iter: " << iter << "\n";
//      std::cout << "eigen: " << eig_new << "\n";
//      std::cout << "delta_e: " << delta_e << "\n";
//      std::cout << "norm: " << norms << "\n";
      eig = eig_new;
      iter++;

    } // end of while loop

    if (iter == max_iter) {
      throw std::runtime_error(
          "Davidson Diagonalization Exceeded Max Iteration");
    }

    return eig;
  };

  /// @return return current eigen vector in Davidson
  value_type &eigen_vector() { return eigen_vector_.back(); }

  // clang-format off
  /**
   *
   * @param HB product with A and guess vector
   * @param B  guess vector
   * @param pred preconditioner, which inherit from DavidsonDiagPred
   *
   * @return B updated guess vector
   * @return updated eigen values, norm of residual
   */
  // clang-format on
  template <typename Pred>
  std::tuple<EigenVector<element_type>, EigenVector<element_type>>
  extrapolate(value_type &HB, value_type &B, const Pred &pred) {
    assert(HB.size() == B.size());

    B_.insert(B_.end(), B.begin(), B.end());
    B.clear();

    HB_.insert(HB_.end(), HB.begin(), HB.end());
    HB.clear();
    // size of new subspace
    const auto n_v = B_.size();

    // compute the new subspace
    compute_new_subspace();

    // do eigen solve
    result_type E(n_roots_);
    RowMatrix<element_type> C(n_v, n_roots_);

    eigen_solve_subspace(E, C);

    // compute eigen_vector at current iteration and store it
    // X(i) = B(i)*C(i)
    value_type X(n_roots_);
    for (std::size_t i = 0; i < n_roots_; ++i) {
      X[i] = copy_and_zero(B_[i]);
      for (std::size_t j = 0; j < n_v; ++j) {
        axpy(X[i], C(j, i), B_[j]);
      }
    }

    // check the size, if exceed n_guess, pop oldest
    if (eigen_vector_.size() == n_guess_) {
      eigen_vector_.pop_front();
    }
    eigen_vector_.push_back(X);

    // compute residual
    // R(i) = (H - e(i)I)*B(i)*C(i)
    //      = (HB(i)*C(i) - e(i)*X(i)
    value_type residual(n_roots_);
    EigenVector<element_type> norms(n_roots_);
    for (std::size_t i = 0; i < n_roots_; ++i) {
      residual[i] = copy_and_zero(X[i]);
      const auto e_i = -E[i];
      axpy(residual[i], e_i, X[i]);
      for (std::size_t j = 0; j < n_v; ++j) {
        axpy(residual[i], C(j, i), HB_[j]);
      }
      norms[i] = norm2(residual[i]);
    }

    // precondition
    // user should define preconditioner
    // usually it is Array(i) = (e(i) - H_D)^-1 R(i)
    // where H_D is the diagonal element of H
    // but H_D can be approximated and computed on the fly
    pred(E, residual);

    // subspace collapse
    // restart with new vector and most recent eigen vector
    // Journal of Computational Chemistry, 11(10), 1164–1168.
    // https://doi.org/10.1002/jcc.540111008
    if (B_.size() > n_roots_ * (max_n_guess_ - 1)) {
      B_.clear();
      HB_.clear();
      subspace_.resize(0, 0);
      B.insert(B.end(), residual.begin(), residual.end());
      // use all stored eigen vector from last n_guess interation
      for (auto &vector : eigen_vector_) {
        B.insert(B.end(), vector.begin(), vector.end());
      }
      // orthognolize all vectors
      gram_schmidt(B);
      // call it second times
      gram_schmidt(B);
    } else {
      // orthognolize new residual with original B
      gram_schmidt(B_, residual);
      // call it twice
      gram_schmidt(B_, residual);
      B = residual;
    }

    return std::make_tuple(E.segment(0, n_roots_), norms);
  }

  /// clean the cached values
  void reset() {
    eigen_vector_.clear();
    HB_.clear();
    B_.clear();
    subspace_.resize(0, 0);
  }

private:

  virtual void compute_new_subspace(){

    // size of original subspace
    const auto n_s = subspace_.cols();
    // size of new subspace
    const auto n_v = B_.size();
    // size of new vector
    const auto n_b = n_v - n_s;

    // compute new subspace
    // G will be replicated Eigen Matrix
    {
      RowMatrix<element_type> G = RowMatrix<element_type>::Zero(n_v, n_v);
      // reuse stored subspace
      G.block(0, 0, n_s, n_s) << subspace_;
      // initialize new value
      for (std::size_t i = 0; i < n_b; ++i) {
        const auto ii = i + n_s;
        for (std::size_t j = 0; j <= ii; ++j) {
          G(ii, j) = dot_product(B_[ii], HB_[j]);
          if (ii != j) {
            G(j, ii) = G(ii, j);
          }
        }
      }

      subspace_ = G;
    }
  }


  virtual void eigen_solve_subspace(result_type& E, RowMatrix<element_type>& C){
    // symmetric matrix
    // this return eigenvalue and eigenvector
    Eigen::SelfAdjointEigenSolver<RowMatrix<element_type>> es(subspace_);

    RowMatrix<element_type> v = es.eigenvectors();
    EigenVector<element_type> e = es.eigenvalues();

    if (es.info() != Eigen::Success) {
      throw std::runtime_error("Eigen::SelfAdjointEigenSolver Failed!\n");
    }

    //        std::cout << es.eigenvalues() << std::endl;

    E = e.segment(0, n_roots_);
    C = v.leftCols(n_roots_);
  }


protected:
  unsigned int n_roots_;
  unsigned int n_guess_;
  unsigned int max_n_guess_;
  std::deque<value_type> eigen_vector_;
  value_type HB_;
  value_type B_;
  RowMatrix<element_type> subspace_;
};

} // namespace code_example
#endif // CODE_EXAMPLE_DAVIDSON_DIAG_H_
