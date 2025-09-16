#include <igl/readOBJ.h>
#include <igl/opengl/glfw/Viewer.h>
#include "Geometry/IO.h"
#include <iostream>

int main() {
    Geometry::GeometryInfo info;
    const std::string root_dir = "C:/Users/Zhehao/Research/SimpleModal/";
    if (!Geometry::IO::readGeometryInfo(root_dir + "asset/plate/plate.json", info)) {
        std::cerr << "Failed to read geometry info from JSON file." << std::endl;
        return -1;
    }

    std::cout << "Geometry Name: " << info.name << std::endl;
    std::cout << "OBJ Path: " << info.obj_path << std::endl;
    std::cout << "Material Name: " << info.material.name << std::endl;
    std::cout << "Material Density: " << info.material.density << std::endl;
    std::cout << "Young's Modulus: " << info.material.youngs_modulus << std::endl;
    std::cout << "Poisson Ratio: " << info.material.poisson_ratio << std::endl;
    std::cout << "Rayleigh Damping Alpha: " << info.material.rayleigh_damping.alpha << std::endl;
    std::cout << "Rayleigh Damping Beta: " << info.material.rayleigh_damping.beta << std::endl;

    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    if (!igl::readOBJ(root_dir + info.obj_path, V, F)) {
        std::cerr << "Failed to load OBJ file." << std::endl;
        return -1;
    }

    igl::opengl::glfw::Viewer viewer;
    viewer.data().set_mesh(V, F);
    viewer.launch();

    return 0;
}
