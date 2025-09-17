#pragma once

#include <vector>
#include <Eigen/Dense>
#include "Geometry/TetMesh.h" // adjust this path as needed
#include <random>

namespace SimpleModal {
class Sampling {

public:
  // Sample N points inside a TetMesh, distributed uniformly by volume
    static void samplePointsInTetMesh(
      Eigen::MatrixXd& points,
      const TetMesh& mesh,
      int num_points,
      unsigned int seed = 42);

    static void samplePointsOnOffsetSurface(
      Eigen::MatrixXd& points,
      const Eigen::MatrixXd& V,
      const Eigen::MatrixXi& F,
      double offset,
      int num_points,
      int offset_surface_resolution,
      unsigned int seed = 1234);
};
}; // namespace SimpleModal
