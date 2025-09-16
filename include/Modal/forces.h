#pragma once

#include <Eigen/Dense>
#include <vector>
#include <tuple>

namespace SimpleModal {

enum class ForceMode {
    CONSTANT,
    GAUSSIAN
};

class ForceTimeSeriesBase
{
public:
    ForceTimeSeriesBase(const double &duration, 
                       const std::vector<Eigen::VectorXd> &force, 
                       const std::vector<std::pair<int, int>> &objectIdx_and_vertexIdx,
                       ForceMode mode = ForceMode::CONSTANT,
                       double gaussian_center_ratio = 0.5,
                       double gaussian_width_ratio = 0.25);

    std::tuple<std::vector<Eigen::VectorXd>, std::vector<std::pair<int, int>>> getForceSampleAtTime(const double &t) const;

    inline int getNumSamples() const { return _num_samples; }

protected:
    const int _sample_rate = 44100;
    int _num_samples;

    double _duration;
    std::vector<std::pair<int, int>> _objectIdx_and_vertexIdx;

    // Force mode and parameters
    ForceMode _force_mode;
    double _gaussian_center_ratio;  // Center of Gaussian pulse as ratio of duration (0-1)
    double _gaussian_radius_ratio;   // Width of Gaussian pulse as ratio of duration (0-1)

    // forces at time sample idx t 
    std::vector<std::vector<Eigen::VectorXd>> _forceBuffer;
};
} // namespace SimpleModal