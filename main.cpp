#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiPlugin.h>
#include <igl/readOBJ.h>
#include <igl/unproject_onto_mesh.h>
#include <igl/per_face_normals.h>

#include "Geometry/FEM.h"
#include "Geometry/Cache.h"
#include "Geometry/EigenSolver.h"
#include "Geometry/IO.h"
#include "Geometry/TetMesh.h"
#include "Geometry/BEM.h"

#include "Transfer/Power.h"
#include "Transfer/EquivalentSources.h"
#include "Transfer/SphericalHarmonics.h"
#include "Transfer/Sampling.h"

#include "Utils/utils.h"
#include "Utils/colormap.h"
#include "Modal/forces.h"
#include "Modal/modalIntegratorIIR.h"
#include "Audio/AudioManager.h"

#include <filesystem>
#include <iostream>
#include <algorithm>
#include "argparse.hpp"

#define _USE_MATH_DEFINES
#include <cmath>

int main(int argc, char *argv[])
{
  std::cout << "Starting main()..." << std::endl;

  // Set up argument parser
  argparse::ArgumentParser program("test_realtime_modal_sound", "1.0");

  // Add command line argument for geometry JSON path
  program.add_argument("geo")
      .help("Name of the geometry")
      .default_value("plate")
      .required();
  program.add_argument("--num_modes", "-m")
      .help("Number of modes to compute")
      .default_value(200)
      .scan<'i', int>();
  program.add_argument("--recompute_modal", "-r")
      .help("Recompute modal data")
      .default_value(false)
      // .default_value(true)
      .implicit_value(true);
  // Parse command line arguments
  try
  {
    program.parse_args(argc, argv);
  }
  catch (const std::runtime_error &err)
  {
    std::cerr << err.what() << std::endl;
    std::cerr << program;
    std::exit(1);
  }

  // Get the geometry JSON path from command line arguments
  std::string geo_name = program.get<std::string>("geo");

  // Load geometry information
  SimpleModal::GeometryInfo info;
  const std::string root_dir = PROJECT_ROOT_DIR "/";
  const std::string asset_dir = root_dir + "asset/";
  if (!SimpleModal::IO::readGeometryInfo(asset_dir + geo_name + "/info.json", info))
  {
    std::cout << "asset_dir: " << asset_dir + geo_name + "/" + geo_name << std::endl;
    std::cerr << "Failed to read geometry info from JSON file." << std::endl;
    return -1;
  }

  // Load the surface mesh
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  if (!igl::readOBJ(root_dir + info.obj_path, V, F))
  {
    std::cerr << "Failed to load OBJ file." << std::endl;
    return -1;
  }
  // Compute per vertex normals
  Eigen::MatrixXd N;
  igl::per_vertex_normals(V, F, N);

  // compute per face normals
  Eigen::MatrixXd per_face_normals;
  igl::per_face_normals(V, F, per_face_normals);

  // Compute per face areas
  Eigen::VectorXd per_face_areas;
  igl::doublearea(V, F, per_face_areas);
  per_face_areas *= 0.5;

  // Load arrow mesh
  Eigen::MatrixXd arrow_V;
  Eigen::MatrixXi arrow_F;
  if (!igl::readOBJ(root_dir + "asset/arrow.obj", arrow_V, arrow_F))
  {
    std::cerr << "Failed to load arrow OBJ file." << std::endl;
    return -1;
  }

  // Center and normalize arrow mesh
  Eigen::Vector3d arrow_center = arrow_V.colwise().mean();
  // arrow_V = arrow_V.rowwise() - arrow_center.transpose();
  double arrow_scale = (arrow_V.rowwise() - arrow_center.transpose()).rowwise().norm().maxCoeff();
  arrow_V /= arrow_scale;

  // Variables for arrow visualization
  bool show_arrow = true;
  Eigen::Vector3d arrow_position;
  Eigen::Vector3d arrow_direction;
  double arrow_size = 0.02; // Adjust this to change arrow size
  // -----------------------------------------------------------------------------------

  // Create cache directory, so we don't need to recompute modal analysis for the same geometry and material every time
  if (!SimpleModal::Cache::createCacheDir(info.obj_path))
  {
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
  int num_modes = program.get<int>("--num_modes");
  bool recompute_modal = program.get<bool>("--recompute_modal");

  bool use_cache = false;

  if (std::filesystem::exists(tet_cache_path) &&
      std::filesystem::exists(modes_cache_path) &&
      !recompute_modal)
  {
    Utils::printInfo("Loading cached tet and mode data from " + tet_cache_path +
                     " and " + modes_cache_path);
    if (SimpleModal::Cache::loadTetMesh(tet_cache_path, tet_mesh) &&
        SimpleModal::Cache::loadModalData(modes_cache_path, S, U))
    {
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
    }
    else
    {
      Utils::printWarning("Failed to load cached data, recomputing...");
    }
  }

  // If we don't use cache, we need to compute the mass and stiffness matrices
  if (!use_cache)
  {
    // Tetrahedralize the mesh
    tet_mesh = SimpleModal::TetMesh::getTetMeshFromSurfaceMesh(root_dir + info.obj_path);

    // Compute mass and stiffness matrices
    TIC(computeMandK);
    Eigen::SparseMatrix<double> M =
        SimpleModal::FEM::computeMassMatrix(tet_mesh, info.material.density);
    Eigen::SparseMatrix<double> K = SimpleModal::FEM::computeStiffnessMatrix(
        tet_mesh, info.material.youngs_modulus, info.material.poisson_ratio);
    TOC(computeMandK);

    // save K and M to be python readable binary files
    std::string K_file = cache_dir + "/" + info.name + "_K.bin";
    std::string M_file = cache_dir + "/" + info.name + "_M.bin";
    SimpleModal::IO::saveMatrixToFile(K_file, K);
    SimpleModal::IO::saveMatrixToFile(M_file, M);

    // Solve the generalized eigenvalue problem
    TIC(EigenSolver);
    num_modes = SimpleModal::EigenSolver::shiftSolve(K, M, U, S, num_modes);
    TOC(EigenSolver);

    // Save the computed data to cache
    if (SimpleModal::Cache::saveTetMesh(tet_cache_path, tet_mesh) &&
        SimpleModal::Cache::saveModalData(modes_cache_path, S, U))
    {
      Utils::printSuccess("Successfully saved data to cache.");
    }
    else
    {
      Utils::printError("Failed to save data to cache.");
    }
  }


  // get surface modal displacement of each mode
  std::vector<Eigen::MatrixXd> u; 
  std::vector<Eigen::MatrixXd> colors;

  for (int mode = 0; mode < num_modes; ++mode)
  {
    Eigen::MatrixXd color = Eigen::MatrixXd::Zero(V.rows(), 3);
    Eigen::VectorXd u_temp = U.col(mode);
    Eigen::MatrixXd full_displacement = u_temp.reshaped(3, tet_mesh.TV.rows());
    Eigen::MatrixXd surface_displacement = Eigen::MatrixXd::Zero(V.rows(), 3);

    for (int i = 0; i < tet_mesh.TV_surface.rows(); ++i)
    {
      int tet_v_idx = tet_mesh.mapSurface2Tet(i);
      surface_displacement.row(i) = full_displacement.col(tet_v_idx);
      color.row(i) = full_displacement.col(tet_v_idx);
    }

    u.push_back(surface_displacement);
    colors.push_back(color);
  }


  // ---------------------------------------------------------------------
  // compute added mass
  // ---------------------------------------------------------------------
  TIC(ComputeAddedMass);
  Eigen::Matrix3d added_mass = Plink::BEM::ComputeAddedMass(V, F, per_face_normals, per_face_areas);
  std::cout << "Added mass: \n"
            << added_mass << std::endl;
  TOC(ComputeAddedMass);





  // ---------------------------------------------------------------------
  // sample points inside mesh 
  // ---------------------------------------------------------------------

  TIC(SamplePoints);

  Eigen::MatrixXd origin = Eigen::MatrixXd::Zero(1, 3);

  int num_sources = 200;
  Eigen::MatrixXd sources(num_sources, 3);
  // SimpleModal::Sampling::samplePointsInTetMesh(sources, tet_mesh, num_sources, 0);
  SimpleModal::Sampling::samplePointsOnOffsetSurface(sources, V, F, 0.3, num_sources, num_sources, 0);

  // Eigen::MatrixXd sources = origin;

  TOC(SamplePoints)

  // ---------------------------------------------------------------------
  



  // ---------------------------------------------------------------------
  // compute transfer coefficients
  // ---------------------------------------------------------------------

  TIC(ComputeTransferCoefficients);

  std::vector<Eigen::VectorXcd> cs =
      Transfer::EquivalentSources::computeTransferCoefficients(V, N, sources, u, S);

  TOC(ComputeTransferCoefficients);
  

  // // ---------------------------------------------------------------------




  // ---------------------------------------------------------------------------------
  // compute equipower dipoles
  // ---------------------------------------------------------------------------------



  TIC(ComputeEquipowerDipoles);

  Eigen::VectorXd powers = Transfer::Power::computePowers(cs, S);
  std::vector<Eigen::Vector3d> equipower_dipole_directions = Transfer::Power::computeEquipowerDipoleDirections(sources, cs, S);

  for (int mode = 0; mode < num_modes; ++mode) {
      std::cout << "-------------------------------------" << std::endl;
      std::cout << "Mode " << mode << std::endl;
      std::cout << "Power: " << powers[mode] << std::endl;
      std::cout << "Equipower Dipole Direction: " << equipower_dipole_directions[mode].transpose() << std::endl;
  }

  TOC(ComputeEquipowerDipoles);


  // ---------------------------------------------------------------------
  // export coefficients for python colormap
  // ---------------------------------------------------------------------


  #include <iomanip>  // for std::setw, std::setprecision

  // make a data/geo_name folder if it doesn't exist
  std::string data_folder = "data/" + geo_name;
  std::filesystem::create_directories(data_folder);

  for(int mode = 0; mode < num_modes; ++mode)
  {
    std::string filename = "data/" + geo_name + "/" + geo_name + "_mode" + std::to_string(mode) + ".csv"; 
    std::ofstream file(filename); 
    if (file.is_open()) 
    {
        file << std::fixed << std::setprecision(6);  // float formatting

        int num_sources = sources.rows(); 
        double k = sqrt(S[mode]) / speed_of_sound; 
        double max_r = V.rowwise().norm().maxCoeff();

        // Write mode and k as the first line
        file << "geo_name," << geo_name << ",mode," << mode << ",k," << k << ",max_r," << max_r << "\n";
        
        // CSV headers
        file << "source_x,source_y,source_z,n,m,re,im\n";

        for (int s = 0; s < num_sources; ++s)
        {
            Eigen::Vector3d source = sources.row(s);
            int offset = 4 * s;

            auto write_row = [&](int n, int m, std::complex<double> coeff) {
                file << source(0) << "," << source(1) << "," << source(2) << ","
                    << n << "," << m << ","
                    << coeff.real() << "," << coeff.imag() << "\n";
            };

            write_row(0,  0, cs[mode][offset + 0]);
            write_row(1, -1, cs[mode](offset + 1));
            write_row(1,  0, cs[mode](offset + 2));
            write_row(1,  1, cs[mode](offset + 3));
        }
        file.close();
    }

    std::string filename_ED = "data/" + geo_name + "/" + geo_name + "_equivalent_dipole_mode" + std::to_string(mode) + ".csv"; 
    std::ofstream file_ED(filename_ED); 
    if (file_ED.is_open()) 
    {
        file_ED << std::fixed << std::setprecision(6);  // float formatting

        double k = sqrt(S[mode]) / speed_of_sound; 
        double max_r = V.rowwise().norm().maxCoeff();

        // Write mode and k as the first line
        file_ED << "geo_name," << geo_name << ",mode," << mode << ",k," << k << ",max_r," << max_r << "\n";

        // CSV headers
        file_ED << "direction_x,direction_y,direction_z,power\n";

        // file_ED << equipower_dipole_directions[mode](0) << "," 
        //         << equipower_dipole_directions[mode](1) << "," 
        //         << equipower_dipole_directions[mode](2) << ","
        //         << powers[mode] << "\n";

        file_ED << equipower_dipole_directions[mode](0) << "," 
                << equipower_dipole_directions[mode](1) << "," 
                << equipower_dipole_directions[mode](2) << ","
                << 1 << "\n";

        file_ED.close();
    }
  }




  // ---------------------------------------------------------------------------------
  // method of known results
  // ---------------------------------------------------------------------------------

  // Utils::printInfo("Method of Known Results");
  // std::cout << "-----------------------------------" << std::endl;

  // std::random_device rd;
  // std::mt19937 gen(rd());
  // // Parameters
  // int num_sample_points = 500;
  // double min_radius = 10.0;

  // // For sample points: range [-50, 50] in each coordinate
  // std::uniform_real_distribution<> coord_dist(-100.0, 100.0);

  // std::vector<Eigen::Vector3d> sample_points;
  // sample_points.reserve(num_sample_points);

  // while (sample_points.size() < num_sample_points) {
  //     Eigen::Vector3d candidate(coord_dist(gen), coord_dist(gen), coord_dist(gen));

  //     // Exclude points inside radius
  //     if (candidate.norm() >= min_radius) {
  //         sample_points.push_back(candidate);
  //     }
  // }

  // for (size_t mode = 0; mode < num_modes; ++mode)
  // {
  //     double k = sqrt(S[mode]) / speed_of_sound;
      
  //     double squared_diff_sum = 0.0;
  //     double squared_exact_sum = 0.0;


  //     int num_sources = sources.rows();

  //     for (const auto &sample_point : sample_points)
  //     {
  //         Eigen::Vector3d sample_source;
  //         sample_source << 0, 0, 0;
  //         // analytical result
  //         double analytical_result = abs(Transfer::SphericalHarmonics::Psi(1, 0, k, sample_point, sample_source));

  //         // compute numerical result
  //         double numerical_result = Transfer::SphericalHarmonics::Pressure(k, sample_point, cs[mode], sources);

  //         // accumulate for RMS
  //         squared_diff_sum  += (numerical_result - analytical_result) * (numerical_result - analytical_result);
  //         squared_exact_sum += analytical_result * analytical_result;

  //         // std::cout << "Relative Error: " << abs(analytical_result - numerical_result) / abs(analytical_result) << std::endl;
  //     }

  //     // compute RMS (L2) error for this mode
  //     double mode_rms = std::sqrt(squared_diff_sum / squared_exact_sum);
  //     std::cout << "Mode " << mode << " RMS error : " << mode_rms << std::endl;
  //     std::cout << "-----------------------------------" << std::endl;
  // }

  // ---------------------------------------------------------------------------------






  // ---------------------------------------------------------------------------------
  // UI 
  // ---------------------------------------------------------------------------------

  // Visualize the eigen modes
  igl::opengl::glfw::Viewer viewer;
  viewer.data().set_mesh(V, F);

  int mode = 0;
  double amplitude = 0.01;
  bool is_picking = true;

  int appliedForceFaceIdx = 0;
  double force_duration = 1.0;
  double force_amplitude = 1.0;

  bool show_wireframe = true;
  bool show_normals = false;

  // Add force mode parameters
  SimpleModal::ForceMode force_mode = SimpleModal::ForceMode::CONSTANT;
  const char *force_mode_names[] = {"Constant", "Gaussian"};
  int current_force_mode = 0;
  double gaussian_center_ratio = 0.5;
  double gaussian_radius_ratio = 0.25;

  auto modalMaterial = SimpleModal::ModalMaterial(info.material.density, info.material.youngs_modulus, info.material.poisson_ratio, info.material.rayleigh_damping.alpha, info.material.rayleigh_damping.beta);
  std::cout << modalMaterial << std::endl;

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

  


  menu.callback_draw_viewer_menu = [&]()
  {
    // Draw parent menu content
    menu.draw_viewer_menu();
    // // Add Wireframe checkbox
    // ImGui::Checkbox("Wireframe", &show_wireframe);
    if (ImGui::CollapsingHeader("Visualize Transfer", ImGuiTreeNodeFlags_DefaultOpen))
    {
      // not working
      if (ImGui::Checkbox("Wireframe", &show_wireframe))
      {
        if(show_wireframe)
        {
          viewer.data().show_faces = false;
        }
        else
        {
          viewer.data().show_faces = true;
        }
      }
      if (ImGui::Checkbox("Show Normals", &show_normals))
      {
        if (show_normals)
        {
            double normal_scale = 0.01; // Adjust based on your mesh size
            Eigen::MatrixXd P1 = V;                         // start points (vertex positions)
            Eigen::MatrixXd P2 = V + normal_scale * N;      // end points (tip of arrows)
            viewer.data(0).add_edges(P1, P2, Eigen::RowVector3d(1, 0, 0)); // red lines
        }
        else
        {
            viewer.data(0).clear_edges(); // <-- This removes the normals
        }
      }
    }
  };
  viewer.data().show_lines = true;
  // viewer.data().show_faces = true;
  viewer.data().show_faces = false;

  // plot source points
  viewer.data(0).add_points(sources, Eigen::RowVector3d(1, 0, 0)); // red points

  // Add these variables at the beginning of main, after other variable declarations
  Utils::Colormap::Type current_colormap = Utils::Colormap::Type::JET;
  const char *colormap_names[] = {"Jet", "UM4", "Bone", "Autumn"};
  int current_colormap_idx = 0;

  
  // Draw additional windows
  menu.callback_draw_custom_window = [&]()
  {
    // Define next window position + size
    ImGui::SetNextWindowPos(ImVec2(180.f * menu.menu_scaling(), 10),
                            ImGuiCond_FirstUseEver);
    ImGui::SetNextWindowSize(ImVec2(300, 400), ImGuiCond_FirstUseEver);
    ImGui::Begin("Simple Modal", nullptr, ImGuiWindowFlags_NoSavedSettings);

    // Display the material name and properties
    ImGui::Text("Material: %s", info.material.name.c_str());
    ImGui::Text("Young's Modulus: %.2f", info.material.youngs_modulus);
    ImGui::Text("Poisson's Ratio: %.2f", info.material.poisson_ratio);
    ImGui::Text("Density: %.2f", info.material.density);
    ImGui::Text("Rayleigh Damping: alpha:%.2f, beta:%.2f", info.material.rayleigh_damping.alpha, info.material.rayleigh_damping.beta);

    // Add colormap selection combo
    if (ImGui::Combo("Colormap", &current_colormap_idx, colormap_names, IM_ARRAYSIZE(colormap_names)))
    {
      current_colormap = static_cast<Utils::Colormap::Type>(current_colormap_idx);
    }

    // Recommand value for semi-implicit Euler: stiffness = 1000, time step =
    // 0.005
    // ImGui::InputInt("Mode", &mode, 1, 100);
    if (ImGui::InputInt("Mode", &mode, 1, 100))
    {
        if (mode < 0) mode = 0;
    }

    ImGui::InputDouble("Amplitude", &amplitude, 0.01, 0.1, "%.4f");


    // Force Parameters Section
    ImGui::Separator();
    ImGui::Checkbox("Is Picking", &is_picking);

    ImGui::Text("Force Parameters");

    // Force mode selection
    if (ImGui::Combo("Force Mode", &current_force_mode, force_mode_names, IM_ARRAYSIZE(force_mode_names)))
    {
      force_mode = static_cast<SimpleModal::ForceMode>(current_force_mode);
    }

    ImGui::InputDouble("Force Duration", &force_duration, 0.01, 0.1, "%.4f");
    ImGui::InputDouble("Force Amplitude", &force_amplitude, 0.01, 0.1, "%.4f");

    // Gaussian parameters (only shown when Gaussian mode is selected)
    if (force_mode == SimpleModal::ForceMode::GAUSSIAN)
    {
      ImGui::PushStyleColor(ImGuiCol_Text, IM_COL32(255, 255, 0, 255)); // Yellow text
      ImGui::Text("Gaussian Parameters");
      ImGui::PopStyleColor();

      ImGui::InputDouble("Center (ratio)", &gaussian_center_ratio, 0.01, 0.1, "%.2f");
      ImGui::SameLine();
      ImGui::TextDisabled("(?)");
      if (ImGui::IsItemHovered())
      {
        ImGui::SetTooltip("Position of Gaussian peak as ratio of duration (0-1)");
      }

      gaussian_center_ratio = std::max(0.0, std::min(1.0, gaussian_center_ratio));

      ImGui::InputDouble("Radius (ratio)", &gaussian_radius_ratio, 0.01, 0.1, "%.2f");
      ImGui::SameLine();
      ImGui::TextDisabled("(?)");
      if (ImGui::IsItemHovered())
      {
        ImGui::SetTooltip("Support Radius of Gaussian pulse as ratio of duration (0-1)");
      }

      gaussian_radius_ratio = std::max(0.01, std::min(1.0, gaussian_radius_ratio));
    }
    ImGui::Text("Arrow Visualization");
    ImGui::Checkbox("Show Arrow", &show_arrow);
    ImGui::InputDouble("Arrow Size", &arrow_size, 0.01, 0.1, "%.3f");

    ImGui::Separator();



    // ------------------------- modal shape visualization ------------------------------
    if (ImGui::Button("Deform", ImVec2(-1, 0)))
    {
      if (mode >= num_modes)
      {
        Utils::printWarning(
            "Mode number is out of range, set to the last mode");
        mode = num_modes - 1;
      }
      Utils::printInfo("Mode " + std::to_string(mode) + " has frequency " +
                       std::to_string(sqrt(S(mode)) / 2 / M_PI) + " Hz");
      // Eigen::VectorXd u = U.col(mode);
      if (amplitude < 1e-10)
      {
        Utils::printWarning("Amplitude is too small, set to 1");
        amplitude = 0.01;
      }

      // Normalize the displacement for color mapping
      Eigen::MatrixXd V_temp = V + u[mode] * amplitude;
      double max_val = colors[mode].rowwise().norm().maxCoeff();
      double min_val = colors[mode].rowwise().norm().minCoeff();

      viewer.data(0).set_vertices(V_temp);
      viewer.data(0).set_colors(Utils::Colormap::mapToColors(
          colors[mode], min_val, max_val, current_colormap));
      viewer.data(0).compute_normals();
    }

    if (ImGui::Button("Reset", ImVec2(-1, 0)))
    {
      viewer.data(0).set_vertices(V);
      viewer.data(0).compute_normals();
    }
    ImGui::End();

    // ImGui::Separator();
  };

  // Initialize audio manager with 44.1kHz sample rate and 1 second buffer
  SimpleModal::AudioManager audioManager(44100, 44100);
  if (!audioManager.initialize())
  {
    std::cerr << "Failed to initialize audio manager" << std::endl;
    return -1;
  }

  // --------------------- apply force with face picking ------------------------------
  viewer.callback_mouse_down =
      [&](igl::opengl::glfw::Viewer &viewer, int, int) -> bool
  {
    if (!is_picking)
    {
      return false;
    }
    int fid;
    Eigen::Vector3f bc;
    // Cast a ray in the view direction starting from the mouse position
    double x = viewer.current_mouse_x;
    double y = viewer.core().viewport(3) - viewer.current_mouse_y;
    if (igl::unproject_onto_mesh(Eigen::Vector2f(x, y), viewer.core().view,
                                 viewer.core().proj, viewer.core().viewport, V, F, fid, bc))
    {
      // paint hit red
      Utils::printInfo("Hit face index = " + std::to_string(fid));
      Eigen::MatrixXd C = Eigen::MatrixXd::Constant(F.rows(), 3, 1);
      C.row(fid) << 1, 0, 0;
      viewer.data(0).set_colors(C);

      // Get face normal and vertex position
      Eigen::MatrixXd FN;
      igl::per_face_normals(V, F, FN);
      Eigen::Vector3d face_normal = FN.row(fid);

      // Get face center for arrow position
      Eigen::Vector3d face_center = (V.row(F(fid, 0)) + V.row(F(fid, 1)) + V.row(F(fid, 2))) / 3.0;

      // Update arrow transformation
      arrow_position = face_center;
      arrow_direction = -face_normal.normalized(); // Opposite to face normal

      // Create rotation matrix to align arrow with direction
      Eigen::Vector3d initial_dir(0, 0, -1);        // Original arrow direction
      Eigen::Vector3d target_dir = arrow_direction; // Direction we want to point to

      // Find rotation axis and angle
      Eigen::Vector3d rotation_axis = initial_dir.cross(target_dir).normalized();
      double cos_angle = initial_dir.dot(target_dir);
      double angle = acos(cos_angle);

      // Create rotation matrix using Rodrigues' formula
      Eigen::Matrix3d rotation;
      if (angle < 1e-6)
      {
        rotation.setIdentity();
      }
      else if (abs(angle - M_PI) < 1e-6)
      {
        // Handle 180-degree rotation case
        Eigen::Vector3d any_perpendicular = initial_dir.cross(Eigen::Vector3d(1, 0, 0));
        if (any_perpendicular.norm() < 1e-6)
        {
          any_perpendicular = initial_dir.cross(Eigen::Vector3d(0, 0, 1));
        }
        rotation_axis = any_perpendicular.normalized();
        rotation = Eigen::AngleAxisd(angle, rotation_axis).matrix();
      }
      else
      {
        rotation = Eigen::AngleAxisd(angle, rotation_axis).matrix();
      }
      // rotation.setIdentity(); // for testing

      // Transform arrow mesh
      Eigen::MatrixXd transformed_arrow = (arrow_V * arrow_size * rotation.transpose()).rowwise() + arrow_position.transpose();
      viewer.data(1).set_vertices(transformed_arrow);
      viewer.data(1).compute_normals();
      viewer.data(1).set_visible(show_arrow);

      // -----------------------------------------------------------------------------
      const int objIdx = 0;
      appliedForceFaceIdx = fid;
      const int appliedForceSurfaceVIdx = F(fid, 0); 
      Utils::printInfo("Hit surface vertex index = " + std::to_string(appliedForceSurfaceVIdx));
      // TODO: apply force to the 3 vertices of the face
      // apply force to the face, force direction is the normal of the face
      Eigen::Vector3d forceDirection = -viewer.data().V_normals.row(appliedForceSurfaceVIdx);

      // transform surface vertex index to tet vertex index
      const int appliedForceTetVIdx = tet_mesh.mapSurface2Tet(appliedForceSurfaceVIdx);
      Eigen::Vector3d force = force_amplitude * forceDirection;

      auto forceRecord = SimpleModal::ForceTimeSeriesBase(
          force_duration,
          {force},
          std::vector<std::pair<int, int>>{std::make_pair(objIdx, appliedForceTetVIdx)},
          force_mode,
          gaussian_center_ratio,
          gaussian_radius_ratio);
      // -----------------------------------------------------------------------------

      const double dt = 1.0 / 44100.0;

      auto modalIntegrator = SimpleModal::ModalIntegratorIIR(dt, S, U, modalMaterial);

      // Prepare audio buffer
      std::vector<float> audioBuffer;
      audioBuffer.reserve(forceRecord.getNumSamples());
      float maxAbsVal = 0.0f;

      TIC(modalIntegration);
      for (int i = 0; i < forceRecord.getNumSamples(); i++)
      {
        auto [forces, objectIdx_and_vertexIdx] = forceRecord.getForceSampleAtTime(i * dt);
        auto audioWaveForm = modalIntegrator.step(i * dt, forces, objectIdx_and_vertexIdx);
        // Store the audio sample
        audioBuffer.push_back(static_cast<float>(audioWaveForm));
        maxAbsVal = std::max(maxAbsVal, std::abs(static_cast<float>(audioWaveForm)));
      }
      TOC(modalIntegration);

      // Set the audio data and start playback
      // Normalize audio buffer to [-1,1] range
      if (maxAbsVal > 0)
      {
        for (float &sample : audioBuffer)
        {
          sample /= maxAbsVal;
        }
      }
      // for (int i = 0; i < audioBuffer.size(); i++)
      // {
      //   std::cout << audioBuffer[i] << ", ";
      // }
      // std::cout << std::endl;
      std::string audioFileName = cache_dir + "/audioBuffer_face" + std::to_string(fid) + ".txt";
      std::ofstream outFile(audioFileName);
      for (int i = 0; i < audioBuffer.size(); i++)
      {
        outFile << audioBuffer[i];
        if (i < audioBuffer.size() - 1)
        {
          outFile << "\n"; // Write each sample on a new line
        }
      }
      outFile.close();
      std::cout << "audioBuffer saved to " << audioFileName << std::endl;

      audioManager.setAudioData(audioBuffer);
      audioManager.resetPlayback();
      audioManager.startPlayback();

      // -----------------------------------------------------------------------------
      auto q_current = modalIntegrator.get_q_current();
      // std::cout << Utils::CYAN << "modal q: " << q_current.transpose() << Utils::RESET << std::endl;

      return true;
    }
    return false;
  };


  // Add arrow mesh as a second object to the viewer
  viewer.append_mesh();
  viewer.data(1).set_mesh(arrow_V, arrow_F);
  viewer.data(1).set_colors(Eigen::RowVector3d(1.0, 0.0, 0.0)); // Red arrow
  viewer.data(1).show_faces = true;
  viewer.data(1).show_lines = false;
  viewer.data(1).set_visible(false); // Initially hidden



  viewer.launch();

  return 0;
}
