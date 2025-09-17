#include "Geometry/EigenSolver.h"
#include "Utils/utils.h"
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/SparseCholesky>

#include <Spectra/GenEigsSolver.h>
#include <Spectra/MatOp/SparseCholesky.h>
#include <Spectra/MatOp/SparseGenMatProd.h>
#include <Spectra/SymGEigsShiftSolver.h>
#include <Spectra/SymGEigsSolver.h>
#include <algorithm>
#include <iostream>

using namespace Spectra;
using namespace SimpleModal;
using namespace std;

class InverseMOperator : public SparseRegularInverse<double>
{
public:
  InverseMOperator(const Eigen::SparseMatrix<double> &M) : SparseRegularInverse<double>(M)
  {
    _diagM = M.diagonal();
  }
  // y = M * x
  void perform_op(const double *x_in, double *y_out) const
  {
    // std::cout << "perform_op" << std::endl;
    Eigen::Map<const Eigen::VectorXd> x(x_in, _diagM.size());
    Eigen::Map<Eigen::VectorXd>(y_out, _diagM.size()) = _diagM.array() * x.array();
  }

  // y = M^{-1} * x
  void solve(const double *x_in, double *y_out) const
  {
    // std::cout << "solve" << std::endl;
    Eigen::Map<const Eigen::VectorXd> x(x_in, _diagM.size());
    Eigen::Map<Eigen::VectorXd>(y_out, _diagM.size()) = x.array() / _diagM.array();
  }

private:
  Eigen::SparseMatrix<double> _M;
  Eigen::VectorXd _diagM;
};

int EigenSolver::solve(const Eigen::SparseMatrix<double> &K,
                       const Eigen::SparseMatrix<double> &M,
                       Eigen::MatrixXd &U, Eigen::VectorXd &S,
                       const int nEigenvalues)
{
  Utils::printInfo("[0] Solving eigen problem using Spectra::SymGEigsSolver");

  // Define matrix operations for K and M
  SparseGenMatProd<double> opK(K);
  // SparseCholesky<double> opM(M);
  // SparseRegularInverse<double> opInvM(M);
  InverseMOperator opInvM(M);

  const int n_dofs = M.rows();
  const int ncv = std::min(2 * nEigenvalues + 1, n_dofs);
  SymGEigsSolver<SparseGenMatProd<double>, InverseMOperator, GEigsMode::RegularInverse>
      solver(opK, opInvM, nEigenvalues, ncv);

  // Initialize and compute
  solver.init();
  int nconv = solver.compute(SortRule::LargestMagn);

  // Retrieve results
  Eigen::VectorXcd evalues;
  if (solver.info() == CompInfo::Successful)
  {
    S = solver.eigenvalues();
    U = solver.eigenvectors();
    std::cout << Utils::CYAN << nconv << " converged eigenvalues found:\n"
              << S.transpose().format(Eigen::IOFormat(
                     Eigen::FullPrecision, 0, ", ", "\n", "", "", "[", "]"))
              << Utils::RESET << std::endl;
  }
  else
  {
    throw std::runtime_error("Eigenvalue decomposition failed");
  }

  return nconv;
}

int EigenSolver::shiftSolve(const Eigen::SparseMatrix<double> &K,
                            const Eigen::SparseMatrix<double> &M,
                            Eigen::MatrixXd &U, Eigen::VectorXd &S,
                            const int nEigenvalues)
{
  Utils::printInfo("[1] Solving eigen problem using Spectra::SymGEigsShiftSolver");
  const int numDofs = M.rows();
  const int femNModes = std::min(nEigenvalues, numDofs - 1);

  /** Compute mass/stiffness eigenvalues and eigenvectors **/
  // if (args.debugMode) {
  // cout << "\nStarting the eigen solver.\n";
  // cout << femNModes << " modes will be computed for FEM analysis.\n\n";
  //}

  // https://spectralib.org/doc/classspectra_1_1symshiftinvert
  // Used to solve y = (K - sigma*M)^{-1}x, which can be optimized
  using OpType = Spectra::SymShiftInvert<double, Eigen::Sparse, Eigen::Sparse>;
  // using BOpType = Spectra::SparseSymMatProd<double>;
  using BOpType = InverseMOperator;

  const int convergenceRatio =
      min(max(2 * femNModes + 1, 20), numDofs);

  const double modesMinFreq = 20.0;    // 20 Hz is the lowest frequency of interest
  const double modesMaxFreq = 20000.0; // 20 kHz is the highest frequency of interest
  const double sigma = pow(2 * M_PI * modesMinFreq, 2);
  const double sigmaMax = pow(2 * M_PI * modesMaxFreq, 2);
  OpType op(K, M);
  BOpType Bop(M);
  Spectra::SymGEigsShiftSolver<OpType, BOpType, Spectra::GEigsMode::ShiftInvert>
      eigs(op, Bop, femNModes, convergenceRatio, sigma);
  eigs.init();
  int nconv = eigs.compute(Spectra::SortRule::LargestMagn, 1000, 1e-10,
                           Spectra::SortRule::SmallestAlge);
  if (eigs.info() == CompInfo::Successful)
  {
    S = eigs.eigenvalues();
    U = eigs.eigenvectors();
    // std::cout << Utils::CYAN << nconv << " converged eigenvalues found:\n"
    //           << S.transpose().format(Eigen::IOFormat(
    //                  Eigen::FullPrecision, 0, ", ", "\n", "", "", "[", "]"))
    //           << utils::RESET << std::endl;

    // TODO: makes this eigenvalue filter a function
    // Filter eigenvalues and eigenvectors within [sigma, sigmaMax]
    std::vector<int> valid_indices;
    for (int i = 0; i < S.size(); ++i)
    {
      if (S[i] >= sigma && S[i] <= sigmaMax)
      {
        valid_indices.push_back(i);
      }
    }

    // Create filtered matrices
    Eigen::VectorXd filtered_S(valid_indices.size());
    Eigen::MatrixXd filtered_U(U.rows(), valid_indices.size());

    for (int i = 0; i < valid_indices.size(); ++i)
    {
      filtered_S[i] = S[valid_indices[i]];
      filtered_U.col(i) = U.col(valid_indices[i]);
    }

    // Replace original matrices with filtered ones
    S = filtered_S;
    U = filtered_U;

    std::cout << Utils::CYAN << "After filtering, " << S.size() 
              << " eigenvalues remain in range [" << sigma << ", " << sigmaMax << "]:\n"
              << S.transpose().format(Eigen::IOFormat(
                     Eigen::FullPrecision, 0, ", ", "\n", "", "", "[", "]"))
              << Utils::RESET << std::endl;

    return valid_indices.size();
  }
  else
  {
    throw std::runtime_error("Eigenvalue decomposition failed");
  }
}
