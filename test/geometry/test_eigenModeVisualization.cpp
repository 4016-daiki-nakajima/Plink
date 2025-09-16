#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiPlugin.h>
#include <igl/readOBJ.h>

#include "FEM.h"
#include "Geometry/EigenSolver.h"
#include "Geometry/IO.h"
#include "Geometry/TetMesh.h"
#include "VegaWrapper.h"
#include "utils.h"
#include <iostream>

#include "VegaWrapper.h"

int main() {
  // Load geometry information
  SimpleModal::GeometryInfo info;
  const std::string root_dir = "/home/zhehaoli/research/SimpleModal/";
  // if (!Geometry::IO::readGeometryInfo(root_dir +
  // "asset/simpleTet/simpleTet.json", info))
  // auto geo_json_path = "asset/bunny/info.json";
  auto geo_json_path = "asset/plate/info.json";
  if (!SimpleModal::IO::readGeometryInfo(root_dir + geo_json_path, info)) {
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

  // Tetrahedralize the mesh
  SimpleModal::TetMesh tet_mesh =
      SimpleModal::TetMesh::getTetMeshFromSurfaceMesh(root_dir + info.obj_path);

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
  Eigen::MatrixXd U;
  Eigen::VectorXd S;
  int num_modes = 15;

  TIC(EigenSolver);
   SimpleModal::EigenSolver::shiftSolve(
   K, M, U, S,
   num_modes); // Solve for the first 3 eigenvalues
  //SimpleModel::EigenSolver::solve(K, M, U, S, num_modes);
  TOC(EigenSolver);

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

  // Draw additional windows
  menu.callback_draw_custom_window = [&]() {
    // Define next window position + size
    ImGui::SetNextWindowPos(ImVec2(180.f * menu.menu_scaling(), 10),
                            ImGuiCond_FirstUseEver);
    ImGui::SetNextWindowSize(ImVec2(200, 160), ImGuiCond_FirstUseEver);
    ImGui::Begin("Eigen Mode", nullptr, ImGuiWindowFlags_NoSavedSettings);

    // Recommand value for semi-implicit Euler: stiffness = 1000, time step =
    // 0.005
    ImGui::InputInt("Mode", &mode, 1, 100);
    ImGui::InputDouble("Amplitude", &amplitude, 0.01, 0.1, "%.4f");

    // ImGui::InputFloat3("(optional) sphere center",
    // mass_spring.sphere_center.data()); Add a button
    if (ImGui::Button("Deform", ImVec2(-1, 0))) {
      if (mode >= num_modes) {
        Utils::printWarning(
            "Mode number is out of range, set to the last mode");
        mode = num_modes - 1;
      }
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
      double max_val = color.cwiseAbs().maxCoeff();
      if (max_val > 0) {
        color /= max_val;
      }
      viewer.data().set_vertices(V_temp);
      viewer.data().set_colors(color);
      viewer.data().compute_normals();
    }
    if (ImGui::Button("Reset", ImVec2(-1, 0))) {
      viewer.data().set_vertices(V);
      viewer.data().compute_normals();
    }
    ImGui::End();
  };

  viewer.launch();

  return 0;
}
