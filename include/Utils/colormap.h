#pragma once
#include <Eigen/Core>
#include <string>

namespace Utils {

class Colormap {
public:
    enum class Type {
        JET,
        UM4,
        BONE,
        AUTUMN
    };

    // Get color for a value in the range [min_value, max_value]
    static Eigen::Vector3f getColor(float value, float min_value, float max_value, Type type = Type::JET);

    // Map a matrix of values to colors
    // Param: `values` is a matrix of size n x 3, each row is a vector of size 3
    static Eigen::MatrixXd mapToColors(const Eigen::MatrixXd& values, double min_value, double max_value, Type type = Type::JET);

private:
    static Eigen::Vector3f interpolateColor(float t, const float colormap[][3], int size);
    static float clamp(float value, float min_value, float max_value);
};

} // namespace SimpleModal
