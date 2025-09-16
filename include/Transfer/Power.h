#pragma once

#include <cmath>
#include <Eigen/Core>
#include "utils.h"



namespace Transfer
{
  class Power
  {
  public:
    static Eigen::VectorXd computePowers(
      const std::vector<Eigen::VectorXcd>& coefficients,
      const Eigen::VectorXd& omega_squareds);

    static std::vector<Eigen::Vector3d> computeEquipowerDipoleDirections(
      const Eigen::MatrixXd &sources,
      const std::vector<Eigen::VectorXcd>& coefficients,
      const Eigen::VectorXd& omega_squareds);

  };
} // namespace Transfer