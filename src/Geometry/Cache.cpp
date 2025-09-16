#include "Geometry/Cache.h"
#include "utils.h"
#include <filesystem>
#include <fstream>
#include <iostream>

namespace SimpleModal {

bool Cache::saveTetMesh(const std::string &filepath,
                        const SimpleModal::TetMesh &tet_mesh) {
  std::ofstream file(filepath, std::ios::binary);
  if (!file.is_open()) {
    Utils::printError("Failed to open file for writing: " + filepath);
    return false;
  }

  // Write number of vertices
  int32_t num_vertices = tet_mesh.TV.rows();
  file.write(reinterpret_cast<const char *>(&num_vertices), sizeof(int32_t));

  // Write vertex positions
  file.write(reinterpret_cast<const char *>(tet_mesh.TV.data()),
             num_vertices * 3 * sizeof(double));

  // Write number of tetrahedrons
  int32_t num_tets = tet_mesh.TT.rows();
  file.write(reinterpret_cast<const char *>(&num_tets), sizeof(int32_t));

  // Write tetrahedron indices
  file.write(reinterpret_cast<const char *>(tet_mesh.TT.data()),
             num_tets * 4 * sizeof(int32_t));
  
  // Write number of surface faces
  int32_t num_faces = tet_mesh.TF.rows();
  file.write(reinterpret_cast<const char *>(&num_faces), sizeof(int32_t));
  file.write(reinterpret_cast<const char *>(tet_mesh.TF.data()),
             num_faces * 3 * sizeof(int32_t));

  return true;
}

bool Cache::loadTetMesh(const std::string &filepath,
                        SimpleModal::TetMesh &tet_mesh) {
  std::ifstream file(filepath, std::ios::binary);
  if (!file.is_open()) {
    Utils::printError("Failed to open file for reading: " + filepath);
    return false;
  }

  // Read number of vertices
  int32_t num_vertices;
  file.read(reinterpret_cast<char *>(&num_vertices), sizeof(int32_t));
  tet_mesh.TV.resize(num_vertices, 3);

  // Read vertex positions
  file.read(reinterpret_cast<char *>(tet_mesh.TV.data()),
            num_vertices * 3 * sizeof(double));

  // Read number of tetrahedrons
  int32_t num_tets;
  file.read(reinterpret_cast<char *>(&num_tets), sizeof(int32_t));
  tet_mesh.TT.resize(num_tets, 4);

  // Read tetrahedron indices
  file.read(reinterpret_cast<char *>(tet_mesh.TT.data()),
            num_tets * 4 * sizeof(int32_t));

  // Read number of surface faces
  int32_t num_faces;
  file.read(reinterpret_cast<char *>(&num_faces), sizeof(int32_t));
  tet_mesh.TF.resize(num_faces, 3);
  file.read(reinterpret_cast<char *>(tet_mesh.TF.data()),
            num_faces * 3 * sizeof(int32_t));
          
  tet_mesh.buildVertexMapping();

  return true;
}

bool Cache::saveModalData(const std::string &filepath, const Eigen::VectorXd &S,
                          const Eigen::MatrixXd &U) {
  std::ofstream file(filepath, std::ios::binary);
  if (!file.is_open()) {
    Utils::printError("Failed to open file for writing: " + filepath);
    return false;
  }

  // Write number of DOFs
  int32_t num_dofs = U.rows();
  file.write(reinterpret_cast<const char *>(&num_dofs), sizeof(int32_t));

  // Write number of modes
  int32_t num_modes = S.size();
  file.write(reinterpret_cast<const char *>(&num_modes), sizeof(int32_t));

  // Write eigenvalues
  file.write(reinterpret_cast<const char *>(S.data()),
             num_modes * sizeof(double));

  // Write eigenvectors
  file.write(reinterpret_cast<const char *>(U.data()),
             num_dofs * num_modes * sizeof(double));

  return true;
}

bool Cache::loadModalData(const std::string &filepath, Eigen::VectorXd &S,
                          Eigen::MatrixXd &U) {
  std::ifstream file(filepath, std::ios::binary);
  if (!file.is_open()) {
    Utils::printError("Failed to open file for reading: " + filepath);
    return false;
  }

  // Read number of DOFs
  int32_t num_dofs;
  file.read(reinterpret_cast<char *>(&num_dofs), sizeof(int32_t));

  // Read number of modes
  int32_t num_modes;
  file.read(reinterpret_cast<char *>(&num_modes), sizeof(int32_t));

  // Read eigenvalues
  S.resize(num_modes);
  file.read(reinterpret_cast<char *>(S.data()), num_modes * sizeof(double));

  // Read eigenvectors
  U.resize(num_dofs, num_modes);
  file.read(reinterpret_cast<char *>(U.data()),
            num_dofs * num_modes * sizeof(double));

  return true;
}

std::string Cache::getCacheDir(const std::string &geometry_name) {
  // Extract the directory path from the geometry name
  size_t last_slash = geometry_name.find_last_of("/");
  std::string base_dir = (last_slash != std::string::npos)
                             ? geometry_name.substr(0, last_slash)
                             : "";

  // Create cache directory path
  return base_dir + "/cache";
}

bool Cache::createCacheDir(const std::string &geometry_name) {
  std::string cache_dir = getCacheDir(geometry_name);

  try {
    if (!std::filesystem::exists(cache_dir)) {
      std::filesystem::create_directories(cache_dir);
      Utils::printInfo("Created cache directory: " + cache_dir);
    } else {
      Utils::printInfo("cache directory exists: " + cache_dir);
    }
    return true;
  } catch (const std::filesystem::filesystem_error &e) {
    std::string error_msg =
        std::string("Failed to create cache directory: ") + e.what();
    Utils::printError(error_msg);
    return false;
  }
}

} // namespace Geometry
