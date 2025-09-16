#include "Geometry/IO.h"
#include "TetMesh.h"
#include "VegaWrapper.h"
#include <Eigen/Dense>
#include <iostream>

int main() {
  Geometry::GeometryInfo info;
  const std::string root_dir = "/home/zhehaoli/research/SimpleModal/";
  if (!Geometry::IO::readGeometryInfo(root_dir + "asset/simpleTet/info.json",
                                      info)) {
    std::cerr << "Failed to read geometry info from JSON file." << std::endl;
    return -1;
  }

  std::cout << "Geometry Name: " << info.name << std::endl;
  std::cout << "OBJ Path: " << info.obj_path << std::endl;
  std::cout << "Material Name: " << info.material.name << std::endl;
  std::cout << "Material Density: " << info.material.density << std::endl;
  std::cout << "Young's Modulus: " << info.material.youngs_modulus << std::endl;
  std::cout << "Poisson Ratio: " << info.material.poisson_ratio << std::endl;
  std::cout << "Rayleigh Damping Alpha: "
            << info.material.rayleigh_damping.alpha << std::endl;
  std::cout << "Rayleigh Damping Beta: " << info.material.rayleigh_damping.beta
            << std::endl;

  // Eigen::MatrixXd V;
  // Eigen::MatrixXi T;
  // if (!igl::readOBJ(root_dir + info.obj_path, V, F)) {
  // std::cerr << "Failed to load OBJ file." << std::endl;
  // return -1;
  //}
  //
  Eigen::MatrixXd TV(4, 3); // Initialize TV as a 4x3 matrix
  TV << 0, 0, 0,            // Set the vertices
      1, 0, 0, 0, 1, 0, 0, 0, 1;

  Eigen::MatrixXi TT(1, 4);
  TT << 0, 1, 2, 3;

  Eigen::MatrixXi TF(4, 3);
  TF << 0, 2, 1, 0, 1, 3, 0, 2, 3, 1, 2, 3;

  SimpleModal::TetMesh tet_mesh(TV, TT, TF);

  // SimpleModal::TetMesh tet_mesh =
  // SimpleModal::TetMesh::getTetMeshFromSurfaceMesh(root_dir + info.obj_path);
  auto [M, K] =
      Vega::StVK::computeMassAndStiffnessMatrices(tet_mesh.TV, tet_mesh.TT, info.material.youngs_modulus, info.material.poisson_ratio, info.material.density);
  std::cout << "M: " << M << std::endl;
  std::cout << "K: " << K << std::endl;
}
