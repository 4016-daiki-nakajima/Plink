#pragma once

#include <Eigen/Eigenvalues>
#include <Eigen/Sparse>
#include "common.h"


namespace SimpleModal
{
  class EigenSolver
  {
  public:
    // Perform generalized eigenvalue decomposition KU = MUS
    static int solve(const Eigen::SparseMatrix<double> &K,
                     const Eigen::SparseMatrix<double> &M, Eigen::MatrixXd &U,
                     Eigen::VectorXd &S, const int nEigenvalues = 3);

    static int shiftSolve(const Eigen::SparseMatrix<double> &K,
                          const Eigen::SparseMatrix<double> &M, Eigen::MatrixXd &U,
                          Eigen::VectorXd &S, const int nEigenvalues = 3);
  };
}; // namespace SimpleModal
