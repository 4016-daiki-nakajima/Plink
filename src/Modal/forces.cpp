#include "forces.h"
#include "utils.h"

using namespace SimpleModal;

ForceTimeSeriesBase::ForceTimeSeriesBase(const double &duration, 
                                       const std::vector<Eigen::VectorXd> &force, 
                                       const std::vector<std::pair<int, int>> &objectIdx_and_vertexIdx,
                                       ForceMode mode,
                                       double gaussian_center_ratio,
                                       double gaussian_radius_ratio)
    : _duration(duration)
    , _objectIdx_and_vertexIdx(objectIdx_and_vertexIdx)
    , _force_mode(mode)
    , _gaussian_center_ratio(gaussian_center_ratio)
    , _gaussian_radius_ratio(gaussian_radius_ratio)
{
    _num_samples = duration * _sample_rate;
    _forceBuffer.resize(_num_samples);

    for (int i = 0; i < _num_samples; i++)
    {
        double t = i * (1.0 / _sample_rate);
        double scale_factor = 1.0;

        if (_force_mode == ForceMode::GAUSSIAN)
        {
            // Calculate Gaussian pulse
            double center_time = _duration * _gaussian_center_ratio;
            double radius = _duration * _gaussian_radius_ratio;
            scale_factor = exp(- (1 / radius) * pow(t - center_time, 2) / pow(_duration, 2));
        }
        // else ForceMode::CONSTANT, scale_factor remains 1.0

        for (int nForce = 0; nForce < force.size(); nForce++)
        {
            _forceBuffer[i].push_back(force[nForce] * scale_factor);
        }
    }
}

std::tuple<std::vector<Eigen::VectorXd>, std::vector<std::pair<int, int>>> ForceTimeSeriesBase::getForceSampleAtTime(const double &t) const
{
    if (t < 0 || t > _duration)
    {
        Utils::printWarning("[ForceTimeSeriesBase] Time is out of range");
    }

    const int sample_idx = (t / _duration) * _num_samples;

    const std::vector<Eigen::VectorXd> forces = _forceBuffer[sample_idx];

    return std::make_tuple(forces, _objectIdx_and_vertexIdx);
}
