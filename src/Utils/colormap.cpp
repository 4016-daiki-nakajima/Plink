#include "Utils/colormap.h"
#include "Utils/colormap_jet.h"
#include "Utils/colormap_bone.h"
#include "Utils/colormap_autumn.h"
#include "Utils/colormap_um4.h"
#include <algorithm>

namespace Utils {

float Colormap::clamp(float value, float min_value, float max_value) {
    return std::min(std::max(value, min_value), max_value);
}

Eigen::Vector3f Colormap::interpolateColor(float t, const float colormap[][3], int size) {
    // Ensure t is in [0,1]
    t = clamp(t, 0.0f, 1.0f);
    
    // Convert t to index
    float index = t * (size - 1);
    int i0 = static_cast<int>(index);
    int i1 = std::min(i0 + 1, size - 1);
    float f = index - i0;  // fractional part

    // Interpolate between the two closest colors
    return Eigen::Vector3f(
        (1 - f) * colormap[i0][0] + f * colormap[i1][0],
        (1 - f) * colormap[i0][1] + f * colormap[i1][1],
        (1 - f) * colormap[i0][2] + f * colormap[i1][2]
    );
}

Eigen::Vector3f Colormap::getColor(float value, float min_value, float max_value, Type type) {
    // Normalize value to [0,1]
    float t = (clamp(value, min_value, max_value) - min_value) / (max_value - min_value);

    switch (type) {
        case Type::JET:
            return interpolateColor(t, colormap_jet, 256);
        case Type::BONE:
            return interpolateColor(t, colormap_bone, 256);
        case Type::AUTUMN:
            return interpolateColor(t, colormap_autumn, 256);
        case Type::UM4:
            return interpolateColor(t, colormap_um4, 1024);
        default:
            return Eigen::Vector3f(t, t, t); // Grayscale as fallback
    }
}

Eigen::MatrixXd Colormap::mapToColors(const Eigen::MatrixXd& values, double min_value, double max_value, Type type) {
    // Create output matrix with same number of rows as input and 3 columns (RGB)
    Eigen::MatrixXd colors(values.rows(), 3);

    // Map each value to a color
    for (int i = 0; i < values.rows(); ++i) {
        // Get color for each value and store in the output matrix
        Eigen::Vector3f color = getColor(values.row(i).norm(), min_value, max_value, type);
        colors.row(i) = color.cast<double>();
    }

    return colors;
}

} // namespace SimpleModal
