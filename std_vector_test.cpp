//
// CopyRight Chong Peng 2017
//

#include <Eigen/Dense>
#include <iostream>

#include "type_traits.h"
#include "std_vector_interface.h"
#include "symm_davidson_diag.h"
#include "nonsymm_davidson_diag.h"

using namespace code_example;

int main() {

  std::cout.precision(12);

  std::cout << "Nonsymmetric Davidson Diagnolization Using std::vector \n";

  // variables
  const std::size_t n = 200;
  const double sparse = 0.1;
  const std::size_t n_roots = 5;   // number of roots to solve in davidson
  const double converge = 1.0e-10; // convergence threshold in davidson
  const std::size_t max_iter = 100; // max iteration in davidson

  // initialize symmetric matrix
  ColMatrix<double> A = ColMatrix<double>::Zero(n, n);
  for (auto i = 0; i < n; i++) {
    A(i, i) = i + 1;
  }
  A = A + sparse * ColMatrix<double>::Random(n, n);

  EigenVector<double> A_diagonal = A.diagonal();
  // eigen solve
  Eigen::EigenSolver<ColMatrix<double>> es(A);

  EigenVector<double> e = es.eigenvalues().real();
  std::sort(e.data(), e.data()+e.size());
  e = e.segment(0,n_roots);

  std::cout << "Reference Result from EigenSolve: " << std::endl
            << e << std::endl;

  /// construct the SymmDavidsonDiag  object
  NonSymmDavidsonDiag<std::vector<double>> dvd(n_roots);

  /// make the initial guess use unit vector
  std::vector<std::vector<double>> guess(n_roots);
  {
    for (std::size_t i = 0; i < n_roots; i++) {
      guess[i] = std::vector<double>(n, 0.0);
      guess[i][i] = 1;
    }
  }

  /// make the preconditioner

  auto pred = [&A_diagonal](const EigenVector<double> &e,
                            std::vector<std::vector<double>> &guess) {

    for (std::size_t i = 0; i < guess.size(); i++) {
      const auto ei = e[i];
      auto &guess_i = guess[i];
      const auto n_r = guess_i.size();
      for (std::size_t i = 0; i < n_r; i++) {
        guess_i[i] = guess_i[i] / (ei - A_diagonal[i]);
      }
    }

  };

  /// make the operator
  auto op = [&A,n](const std::vector<std::vector<double>> &vec) {
    const std::size_t n_vec = vec.size();

    std::vector<std::vector<double>> HC(n_vec);

    const char trans = 'N';
    const int32_t rows = A.rows();
    const int32_t cols = A.cols();
    const double alpha = 1.0;
    const double beta = 0.0;
    const int32_t inc = 1;
    for (std::size_t i = 0; i < n_vec; i++) {
      HC[i] = std::vector<double>(vec[i].size(), 0.0);
      dgemv_(&trans, &rows, &cols, &alpha, A.data(), &rows, vec[i].data(), &inc,
             &beta, HC[i].data(), &inc);
    }

    return HC;
  };

  /// solve

  auto eig = dvd.solve(guess, op, pred, converge, max_iter);

  std::cout << "NonSymmDavidsonDiag result: " << std::endl;
  std::cout << eig << std::endl;

  return 0;
}
