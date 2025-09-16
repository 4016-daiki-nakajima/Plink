#pragma once

#include <Eigen/Dense>
#include <vector>
#include "modalData.h"

namespace SimpleModal {

class ModalIntegratorBase {
public:
  ModalIntegratorBase(const double dt, const Eigen::VectorXd &eigenValues,
                      const Eigen::MatrixXd &eigenVectors, 
                      const ModalMaterial &material)
      : _dt(dt), _eigenValues(eigenValues), _eigenVectors(eigenVectors),
        _material(material)    
  {};
  virtual double step(double time,
                      const std::vector<Eigen::VectorXd> &forces,
                      const std::vector<std::pair<int, int>> &objectIdx_and_vertexIdx) = 0;

protected:
  double _dt = 0;
  Eigen::VectorXd _eigenValues;
  Eigen::MatrixXd _eigenVectors;
  ModalMaterial _material;
};

}; // namespace SimpleModal
