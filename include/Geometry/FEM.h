#pragma once

#include "TetMesh.h"
#include <Eigen/Sparse>
#include "common.h"

namespace SimpleModal {

class FEM {
public:

  double static computeTetrahedronVolume(const Eigen::MatrixXd &TV, const Eigen::Vector4i &tet);
  // Compute the lumped mass matrix for a given tetrahedral mesh
  static Eigen::SparseMatrix<double>
  computeMassMatrix(const TetMesh &mesh, const double &density = 1.0);

  // Compute the stiffness matrix for a given tetrahedral mesh
  static Eigen::SparseMatrix<double>
  computeStiffnessMatrix(const TetMesh &mesh, const double &young_modulus = 1E7,
                         const double &poisson_ratio = 0.3);
};
}; // namespace SimpleModal
