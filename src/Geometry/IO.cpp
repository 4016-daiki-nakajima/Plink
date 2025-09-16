#include "Geometry/IO.h"
#include "Utils/utils.h"
#include <fstream>
#include "json.hpp" // Include the JSON library

namespace SimpleModal
{

    bool IO::readGeometryInfo(const std::string &json_file, GeometryInfo &info)
    {
        std::ifstream file(json_file);
        if (!file.is_open())
        {
            return false;
        }

        nlohmann::json json_data;
        file >> json_data;

        info.name = json_data["name"];
        info.obj_path = json_data["obj_path"];
        info.material.name = json_data["material"]["name"];
        info.material.density = json_data["material"]["density"];
        info.material.youngs_modulus = json_data["material"]["youngs_modulus"];
        info.material.poisson_ratio = json_data["material"]["poisson_ratio"];
        info.material.rayleigh_damping.alpha = json_data["material"]["Rayleigh_damping"]["alpha"];
        info.material.rayleigh_damping.beta = json_data["material"]["Rayleigh_damping"]["beta"];
        try
        {
            info.modal.assumed_num_modes = json_data["Modal"]["num_modes"];
            info.modal.assumed_max_freq = json_data["Modal"]["max_freq"];
        }
        catch (const std::exception &e)
        {
            Utils::printWarning("No modal info found in " + json_file + ". Using default values.");
            info.modal.assumed_num_modes = 100;
            info.modal.assumed_max_freq = 20000.0;
        }

        return true;
    }

    // Function to save a sparse matrix to a binary file in a Python-readable format
    bool IO::saveMatrixToFile(const std::string &filename, const Eigen::SparseMatrix<double> &matrix)
    {
        std::ofstream file(filename, std::ios::binary);
        if (!file.is_open())
        {
            Utils::printError("Failed to open file for writing: " + filename);
            return false;
        }

        // Write matrix dimensions
        int rows = matrix.rows();
        int cols = matrix.cols();
        int nnz = matrix.nonZeros();
        file.write(reinterpret_cast<const char *>(&rows), sizeof(int));
        file.write(reinterpret_cast<const char *>(&cols), sizeof(int));
        file.write(reinterpret_cast<const char *>(&nnz), sizeof(int));

        // Write the outer index pointers (column starts)
        file.write(reinterpret_cast<const char *>(matrix.outerIndexPtr()), (cols + 1) * sizeof(int));

        // Write the inner indices (row indices)
        file.write(reinterpret_cast<const char *>(matrix.innerIndexPtr()), nnz * sizeof(int));

        // Write the values
        file.write(reinterpret_cast<const char *>(matrix.valuePtr()), nnz * sizeof(double));

        return true;
    }

} // namespace SimpleModal
