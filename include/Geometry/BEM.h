#pragma once

#include <Eigen/Core>

namespace Plink
{
  class BEM
  {
  public:
    static Eigen::Matrix3d ComputeAddedMass(
        const Eigen::MatrixXd &V,
        const Eigen::MatrixXi &F,
        const Eigen::MatrixXd &per_face_normals,
        const Eigen::VectorXd &per_face_areas);

  };
} // namespace Plink
