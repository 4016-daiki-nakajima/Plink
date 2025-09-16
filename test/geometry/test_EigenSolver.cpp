#include "EigenSolver.h"
#include "FEM.h"
#include "TetMesh.h"
#include "utils.h"
#include <iostream>

using namespace SimpleModal;

int main() {
  // Load a tetrahedral mesh
  const std::string root_path = "C:/Users/Zhehao/Research/SimpleModal/";
  TetMesh mesh =
      TetMesh::getTetMeshFromSurfaceMesh(root_path + "asset/fertility.off");
  // TetMesh mesh = TetMesh::getTetrahedralMesh("asset/stanford-bunny.obj");

  // Compute mass and stiffness matrices
  TIC(computeMassMatrix);
  Eigen::SparseMatrix<double> M = FEM::computeMassMatrix(mesh);
  TOC(computeMassMatrix);

  TIC(computeStiffnessMatrix);
  Eigen::SparseMatrix<double> K = FEM::computeStiffnessMatrix(mesh);
  TOC(computeStiffnessMatrix);

  // Perform generalized eigenvalue decomposition
  Eigen::MatrixXd U;
  Eigen::VectorXd S;
  try {
    TIC(eigenSolver);
    EigenSolver::solve(K, M, U, S, 3);
    TOC(eigenSolver);
  } catch (const std::exception &e) {
    std::cerr << "Error: " << e.what() << std::endl;
    return -1;
  }

  // Verify the results
  Eigen::MatrixXd I = U.transpose() * M * U;
  Eigen::MatrixXd S_check = U.transpose() * K * U;

  std::cout << std::endl
            << Utils::GREEN << "U.transpose() * M * U (should be identity):\n"
            << I << Utils::RESET << std::endl;
  std::cout << std::endl
            << Utils::GREEN << "U.transpose() * K * U (should be diagonal):\n"
            << S_check << Utils::RESET << std::endl;

  return 0;
}
