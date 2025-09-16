#include "Geometry/FEM.h"
#include "Geometry/TetMesh.h"
#include "utils.h"
#include <Eigen/Dense>

using namespace SimpleModal;

int main() {
  Eigen::MatrixXd TV(4, 3); // Initialize TV as a 4x3 matrix
  TV << 0, 0, 0,            // Set the vertices
      1, 0, 0, 0, 1, 0, 0, 0, 1;

  Eigen::MatrixXi TT(1, 4);
  TT << 0, 1, 2, 3;

  Eigen::MatrixXi TF(4, 3);
  TF << 0, 2, 1, 0, 1, 3, 0, 2, 3, 1, 2, 3;

  TetMesh tet_mesh(TV, TT, TF);

  // Compute the mass matrix
  Eigen::SparseMatrix<double> M = FEM::computeMassMatrix(tet_mesh, 1.0);

  // Compute the stiffness matrix
  Eigen::SparseMatrix<double> K =
      FEM::computeStiffnessMatrix(tet_mesh, 1, 0.45);

  std::cout << "M: " << M << std::endl;

  std::cout << "K: " << K << std::endl;

  return 0;
}
