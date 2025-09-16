#pragma once

#include <vector>
#include <cmath>
#include "utils.h"
namespace SimpleModal
{
  class ModaData
  {
  public:
    ModaData() {};
    ModaData(std::vector<double> omegaSquared,
             std::vector<std::vector<double>> modes)
        : omegaSquared(omegaSquared), modes(modes) {};

    std::vector<double> omegaSquared;
    std::vector<std::vector<double>> modes;

    int N_modesAudible = -1;
    double freqThresCache = 22100;

    // add function from IO.h here to load modal data from a file
  };

  class ModalMaterial
  {
  public:
    ModalMaterial() {};
    ModalMaterial(double density, double youngsModulus, double poissonRatio, double alpha, double beta)
        : density(density), youngsModulus(youngsModulus), poissonRatio(poissonRatio), alpha(alpha), beta(beta) {};

    friend std::ostream &operator<<(std::ostream &os, const ModalMaterial &material)
    {
      os << Utils::GREEN << "[ModalMaterial] density: " << material.density << ", youngsModulus: " << material.youngsModulus << ", poissonRatio: " << material.poissonRatio << ", alpha: " << material.alpha << ", beta: " << material.beta << Utils::RESET;
      return os;
    }

    double density = 0;
    double alpha = 0;
    double beta = 0;
    double youngsModulus = 0;
    double poissonRatio = 0;

    // formula reference: https://www.cs.cornell.edu/~djames/papers/DyRT.pdf
    inline double one_minus_nu2_over_E() const
    {
      return (1.0 - std::pow(poissonRatio, 2)) / youngsModulus;
    }
    inline double xi(const double &omega_i) const
    {
      return 0.5 * (alpha / omega_i + beta * omega_i);
    } // DyRT eq. 10, xi = [0, 1]
    inline double omega_di(const double &omega_i) const
    {
      return omega_i * sqrt(1.0 - pow(xi(omega_i), 2));
    } // DyRT eq. 12.
  };
}; // namespace SimpleModal
