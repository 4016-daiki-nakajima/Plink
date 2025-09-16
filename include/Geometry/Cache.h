#pragma once

#include <string>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include "TetMesh.h"
#include "common.h"

namespace SimpleModal
{

    class Cache
    {
    public:
        // Save tetrahedral mesh to a binary file
        static bool saveTetMesh(const std::string &filepath, const SimpleModal::TetMesh &tet_mesh);

        // Load tetrahedral mesh from a binary file
        static bool loadTetMesh(const std::string &filepath, SimpleModal::TetMesh &tet_mesh);

        // Save modal data (eigenvalues and eigenvectors) to a binary file
        static bool saveModalData(const std::string &filepath,
                                  const Eigen::VectorXd &S,
                                  const Eigen::MatrixXd &U);

        // Load modal data from a binary file
        static bool loadModalData(const std::string &filepath,
                                  Eigen::VectorXd &S,
                                  Eigen::MatrixXd &U);

        // Get the cache directory path for a given geometry name
        static std::string getCacheDir(const std::string &geometry_name);

        // Create cache directory if it doesn't exist
        static bool createCacheDir(const std::string &geometry_name);
    };

} // namespace Geometry