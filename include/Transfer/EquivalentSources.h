#pragma once

#include <cmath>
#include <Eigen/Core>


namespace Transfer
{
  class EquivalentSources
  {
  public:

    static std::vector<Eigen::VectorXcd> computeTransferCoefficients(
        const Eigen::MatrixXd &V,
        const Eigen::MatrixXd &N,
        const Eigen::MatrixXd &sources,
        const std::vector<Eigen::MatrixXd> &displacements,
        Eigen::VectorXd omega_squareds);
  };
} // namespace Transfer