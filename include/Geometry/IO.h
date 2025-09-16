#pragma once

#include <string>
#include <Eigen/Core>
#include <Eigen/Sparse>

namespace SimpleModal {

struct RayleighDamping {
    double alpha;
    double beta;
};

struct Material {
    std::string name;
    double density;
    double youngs_modulus;
    double poisson_ratio;
    RayleighDamping rayleigh_damping;
};

struct ModalInfo {
    int assumed_num_modes;
    double assumed_max_freq;
};

struct GeometryInfo {
    std::string name;
    std::string obj_path;
    Material material;
    ModalInfo modal;
};

class IO {
public:
    static bool readGeometryInfo(const std::string& json_file, GeometryInfo& info);
    static bool saveMatrixToFile(const std::string& filename, const Eigen::SparseMatrix<double>& matrix);

};

} // namespace SimpleModal
