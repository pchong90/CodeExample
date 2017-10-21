//
// CopyRight Chong Peng 2017
//

#include <Eigen/Dense>
#include <iostream>

#include "eigen_interface.h"
#include "type_traits.h"
#include "symm_davidson_diag.h"
#include "nonsymm_davidson_diag.h"

using namespace code_example;

int main() {

  std::cout.precision(12);

  std::cout << "Davidson Diagnolization Using Eigen Matrix \n";
  
  /// typedef the Array type to use
  using Array = RowMatrix<double>;

  // variables
  const std::size_t n = 500;
  const double sparse = 0.1;
  const std::size_t n_roots = 5;   // number of roots to solve in davidson
  const double converge = 1.0e-10; // convergence threshold in davidson
  const std::size_t max_iter = 40; // max iteration in davidson

  // initialize symmetric matrix
  RowMatrix<double> A = RowMatrix<double>::Zero(n, n);
  for (auto i = 0; i < n; i++) {
    A(i, i) = i + 1;
  }
  A = A + sparse * RowMatrix<double>::Random(n, n);

  RowMatrix<double> A_T = A.transpose();
  A = 0.5 * (A_T + A);

  EigenVector<double> A_diagonal = A.diagonal();
  // eigen solve
  Eigen::SelfAdjointEigenSolver<RowMatrix<double>> es(A);
  EigenVector<double> e = es.eigenvalues().segment(0, n_roots);

  std::cout << "EigenSolve result: " << std::endl << e << std::endl;

  /// construct the SymmDavidsonDiag  object
  SymmDavidsonDiag<Array> dvd(n_roots);

  /// make the initial guess use unit vector
  std::vector<Array> guess(n_roots);
  {
    Array guess_all = Array::Identity(n, n_roots);
    for(std::size_t i = 0; i < n_roots; i++){
      guess[i] = guess_all.col(i);
    }
  }

  /// make the preconditioner

  auto pred = [&A_diagonal](const EigenVector<double> &e,
                            std::vector<Array> &guess) {

    for (std::size_t i = 0; i < guess.size(); i++) {
      const auto ei = e[i];
      auto& guess_i = guess[i];
      assert(guess_i.cols() == 1);
      const auto n_r = guess_i.rows();
      for(std::size_t i = 0; i < n_r; i++){
        guess_i(i,0) = guess_i(i,0)/ (ei - A_diagonal[i]);
      }
    }

  };

  /// make the operator
  auto op = [& A](const std::vector<Array>& vec){
    const std::size_t n = vec.size();

    std::vector<Array> HC(n);

    for(std::size_t i = 0; i < n; i++){
      HC[i] = A*vec[i];
    }

    return HC;
  };

  /// solve

  auto eig = dvd.solve(guess, op, pred, converge, max_iter);

  std::cout << "SymmDavidsonDiag result: " << std::endl;
  std::cout << eig << std::endl;

  return 0;
}
