#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiPlugin.h>
#include <igl/readOBJ.h>

#include "FEM.h"
#include "Geometry/Cache.h"
#include "Geometry/EigenSolver.h"
#include "Geometry/IO.h"
#include "Geometry/TetMesh.h"
#include "utils.h"
#include <filesystem>
#include <iostream>

#include "argparse.hpp"

int main(int argc, char* argv[]) {
  // Set up argument parser
  argparse::ArgumentParser program("test_cache", "1.0");
  
  // Add command line argument for geometry JSON path
  program.add_argument("--geo", "-g")
    .help("Name of the geometry")
    .default_value("plate")
    .required();

  // Parse command line arguments
  try {
    program.parse_args(argc, argv);
  } catch (const std::runtime_error& err) {
    std::cerr << err.what() << std::endl;
    std::cerr << program;
    std::exit(1);
  }

  // Get the geometry JSON path from command line arguments
  std::string geo_name = program.get<std::string>("--geo");

  // Load geometry information
  Geometry::GeometryInfo info;
  const std::string root_dir = "/home/zhehaoli/research/SimpleModal/";
  if (!Geometry::IO::readGeometryInfo(root_dir + "asset/" + geo_name + "/info.json", info)) {
    std::cerr << "Failed to read geometry info from JSON file." << std::endl;
    return -1;
  }

  // Load the surface mesh
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  if (!igl::readOBJ(root_dir + info.obj_path, V, F)) {
    std::cerr << "Failed to load OBJ file." << std::endl;
    return -1;
  }

  // Create cache directory
  if (!Geometry::Cache::createCacheDir(info.obj_path)) {
    std::cerr << "Failed to create cache directory." << std::endl;
    return -1;
  }

  // Try to load cached data
  std::string cache_dir = Geometry::Cache::getCacheDir(info.obj_path);
  std::string tet_cache_path = cache_dir + "/tet.bin";
  std::string modes_cache_path = cache_dir + "/modes.bin";

  SimpleModal::TetMesh tet_mesh;
  Eigen::MatrixXd U;
  Eigen::VectorXd S;
  int num_modes = 15;

  bool use_cache = false;

  if (std::filesystem::exists(tet_cache_path) &&
      std::filesystem::exists(modes_cache_path)) {
    Utils::printInfo("Loading cached tet and mode data from " + tet_cache_path +
                     " and " + modes_cache_path);
    if (Geometry::Cache::loadTetMesh(tet_cache_path, tet_mesh) &&
        Geometry::Cache::loadModalData(modes_cache_path, S, U)) {
      use_cache = true;
      Utils::printSuccess("Successfully loaded cached data.");
      Utils::printInfo("Number of tet vertices: " +
                       std::to_string(tet_mesh.TV.rows()));
      Utils::printInfo("Number of tet tets: " +
                       std::to_string(tet_mesh.TT.rows()));
      Utils::printInfo("Number of tet surface faces: " +
                       std::to_string(tet_mesh.TF.rows()));
      Utils::printInfo("Number of modes: " + std::to_string(S.size()));
    } else {
      Utils::printWarning("Failed to load cached data, recomputing...");
    }
  }

  if (!use_cache) {
    // Tetrahedralize the mesh
    tet_mesh = SimpleModal::TetMesh::getTetMeshFromSurfaceMesh(root_dir +
                                                               info.obj_path);

    // Compute mass and stiffness matrices
    TIC(computeMandK);
    Eigen::SparseMatrix<double> M =
        SimpleModel::FEM::computeMassMatrix(tet_mesh, info.material.density);
    Eigen::SparseMatrix<double> K = SimpleModel::FEM::computeStiffnessMatrix(
        tet_mesh, info.material.youngs_modulus, info.material.poisson_ratio);

    // auto [M, K] = Vega::StVK::computeMassAndStiffnessMatrices(
    // tet_mesh.TV, tet_mesh.TT, info.material.youngs_modulus,
    // info.material.poisson_ratio, info.material.density);
    TOC(computeMandK);

    // Solve the generalized eigenvalue problem
    TIC(EigenSolver);
    SimpleModal::EigenSolver::shiftSolve(K, M, U, S, num_modes);
    TOC(EigenSolver);

    // Save the computed data to cache
    if (Geometry::Cache::saveTetMesh(tet_cache_path, tet_mesh) &&
        Geometry::Cache::saveModalData(modes_cache_path, S, U)) {
      Utils::printSuccess("Successfully saved data to cache.");
    } else {
      Utils::printError("Failed to save data to cache.");
    }
  }

  // Visualize the eigen modes
  igl::opengl::glfw::Viewer viewer;
  viewer.data().set_mesh(V, F);

  int mode = 0;
  double amplitude = 0.01;