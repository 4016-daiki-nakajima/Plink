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
#include "Utils/colormap.h"
#include <filesystem>
#include <iostream>
// #include <cmath>
#define _USE_MATH_DEFINES
#include <cmath>


#include "argparse.hpp"

int main(int argc, char* argv[]) {
  // Set up argument parser
  argparse::ArgumentParser program("test_cache", "1.0");
  
  // Add command line argument for geometry JSON path
  program.add_argument("--geo", "-g")
    .help("Name of the geometry")
    .default_value("plate")
    .required();

  program.add_argument("--no-cache", "-nc")
    .help("Skip using and storing cache")
    .default_value(false)
    .implicit_value(true);

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
  bool user_explicit_no_cache = program.get<bool>("--no-cache");

  // Load geometry information
  SimpleModal::GeometryInfo info;
  const std::string root_dir = "/home/zhehaoli/research/SimpleModal/";
  if (!SimpleModal::IO::readGeometryInfo(root_dir + "asset/" + geo_name + "/info.json", info)) {
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
  if (!SimpleModal::Cache::createCacheDir(info.obj_path)) {
    std::cerr << "Failed to create cache directory." << std::endl;
    return -1;
  }

  // Try to load cached data
  std::string cache_dir = SimpleModal::Cache::getCacheDir(info.obj_path);
  std::string tet_cache_path = cache_dir + "/tet.bin";
  std::string modes_cache_path = cache_dir + "/modes.bin";

  SimpleModal::TetMesh tet_mesh;
  Eigen::MatrixXd U;
  Eigen::VectorXd S;
  int num_modes = std::max(info.modal.assumed_num_modes+10, 100);
  Utils::printInfo("Assuming " + std::to_string(num_modes) + " modes.");

  bool use_cache = false;

  if (!user_explicit_no_cache && std::filesystem::exists(tet_cache_path) &&
      std::filesystem::exists(modes_cache_path)) {
    Utils::printInfo("Loading cached tet and mode data from " + tet_cache_path +
                     " and " + modes_cache_path);
    if (SimpleModal::Cache::loadTetMesh(tet_cache_path, tet_mesh) &&
        SimpleModal::Cache::loadModalData(modes_cache_path, S, U)) {
      use_cache = true;
      num_modes = S.size();
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
        SimpleModal::FEM::computeMassMatrix(tet_mesh, info.material.density);
    Eigen::SparseMatrix<double> K = SimpleModal::FEM::computeStiffnessMatrix(
        tet_mesh, info.material.youngs_modulus, info.material.poisson_ratio);

    // auto [M, K] = Vega::StVK::computeMassAndStiffnessMatrices(
    // tet_mesh.TV, tet_mesh.TT, info.material.youngs_modulus,
    // info.material.poisson_ratio, info.material.density);
    TOC(computeMandK);

    // Solve the generalized eigenvalue problem
    TIC(EigenSolver);
    num_modes = SimpleModal::EigenSolver::shiftSolve(K, M, U, S, num_modes);
    // SimpleModel::EigenSolver::solve(K, M, U, S, num_modes);
    TOC(EigenSolver);

    // Save the computed data to cache
    // if (!user_explicit_no_cache) {
      if (SimpleModal::Cache::saveTetMesh(tet_cache_path, tet_mesh) &&
          SimpleModal::Cache::saveModalData(modes_cache_path, S, U)) {
        Utils::printSuccess("Successfully saved data to cache.");
      } else {
        Utils::printError("Failed to save data to cache.");
      }
    // }
  }

  // Visualize the eigen modes
  igl::opengl::glfw::Viewer viewer;
  viewer.data().set_mesh(V, F);

  int mode = 0;
  double amplitude = 0.01;

  // ------------------------------ basic setting ------------------------------
  // Init the viewer
  // Attach a menu plugin
  igl::opengl::glfw::imgui::ImGuiPlugin plugin;
  // igl::opengl::glfw::imgui::ImGuizmoPlugin guizmo_plugin;
  igl::opengl::glfw::imgui::ImGuiMenu menu;
  viewer.plugins.push_back(&plugin);
  // viewer.plugins.push_back(&guizmo_plugin);
  plugin.widgets.push_back(&menu);

  // Add content to the default menu window
  menu.callback_draw_viewer_menu = [&]() {
    // Draw parent menu content
    menu.draw_viewer_menu();
  };

  // Add these variables at the beginning of main, after other variable declarations
  Utils::Colormap::Type current_colormap = Utils::Colormap::Type::JET;
  const char* colormap_names[] = {"JET", "UM4", "BONE", "AUTUMN"};
  int current_colormap_idx = 0;

  // Draw additional windows
  menu.callback_draw_custom_window = [&]() {
    // Define next window position + size
    ImGui::SetNextWindowPos(ImVec2(180.f * menu.menu_scaling(), 10),
                            ImGuiCond_FirstUseEver);
    ImGui::SetNextWindowSize(ImVec2(200, 200), ImGuiCond_FirstUseEver);
    ImGui::Begin("Eigen Mode", nullptr, ImGuiWindowFlags_NoSavedSettings);

    // Add colormap selection combo
    if (ImGui::Combo("Colormap", &current_colormap_idx, colormap_names, IM_ARRAYSIZE(colormap_names))) {
        current_colormap = static_cast<Utils::Colormap::Type>(current_colormap_idx);
    }

    // Recommand value for semi-implicit Euler: stiffness = 1000, time step = 0.005
    ImGui::InputInt("Mode", &mode, 1, 100);
    ImGui::InputDouble("Amplitude", &amplitude, 0.01, 0.1, "%.4f");

    if (ImGui::Button("Deform", ImVec2(-1, 0))) {
      if (mode >= num_modes) {
        Utils::printWarning(
            "Mode number is out of range, set to the last mode");
        mode = num_modes - 1;
      }
      Utils::printInfo("Mode " + std::to_string(mode) + " has frequency " +
                       std::to_string(sqrt(S(mode)) / 2 / M_PI) + " Hz");
      Eigen::VectorXd u = U.col(mode);
      // this displacement contains inner vertices of the tetrahedral mesh
      Eigen::MatrixXd displacement =
          u.reshaped(3, tet_mesh.TV.rows()) * amplitude;
      Eigen::MatrixXd V_temp = V;
      std::cout << Utils::YELLOW
                << "max displacement: " << displacement.cwiseAbs().maxCoeff()
                << Utils::RESET << std::endl;

      Eigen::MatrixXd color = Eigen::MatrixXd::Zero(V.rows(), 3);
      // traverse all surface triangles
      for (int f = 0; f < tet_mesh.TF.rows(); f++) {
        for (int i = 0; i < 3; i++) {
          auto v_idx = tet_mesh.TF(f, i);
          if (V_temp.row(v_idx) == V.row(v_idx)) // not deformed yet
          {
            V_temp.row(v_idx) += displacement.col(v_idx);
            color.row(v_idx) = displacement.col(v_idx);
          }
        }
      }
      // Normalize the displacement for color mapping
      double max_val = color.rowwise().norm().maxCoeff();
      double min_val = color.rowwise().norm().minCoeff();
      std::cout << "color max_val: " << max_val << std::endl;
      viewer.data().set_vertices(V_temp);
      // Update color mapping with selected colormap
      viewer.data().set_colors(Utils::Colormap::mapToColors(
          color, min_val, max_val, current_colormap));
      viewer.data().compute_normals();
    }
    if (ImGui::Button("Reset", ImVec2(-1, 0))) {
      viewer.data().set_vertices(V);
      viewer.data().set_colors(Eigen::MatrixXd::Ones(V.rows(), 3));  // Reset to white
      viewer.data().compute_normals();
    }
    ImGui::End();
  };

  viewer.launch();

  return 0;
}
