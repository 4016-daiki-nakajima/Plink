#include "Modal/modalIntegratorIIR.h"

using namespace SimpleModal;

ModalIntegratorIIR::ModalIntegratorIIR(const double dt,
                                       const Eigen::VectorXd &eigenValues,
                                       const Eigen::MatrixXd &eigenVectors,
                                       const ModalMaterial &material)
    : ModalIntegratorBase(dt, eigenValues, eigenVectors, material)
{

  _N_modes = eigenValues.size();
  _q_prev.setZero(_N_modes);
  _q_curr.setZero(_N_modes);
  _q_next.setZero(_N_modes);
  _q_next_next.setZero(_N_modes);

  // Compute filter coefficients for time step
  _coeff_qNew.setZero(_N_modes);
  _coeff_qOld.setZero(_N_modes);
  _coeff_Q.setZero(_N_modes);

  for (int i = 0; i < _N_modes; i++)
  {
    // FIXME: What if _eigenValues[i] < 0? 
    if (_eigenValues[i] < 0) {
      Utils::printError("[ModalIntegratorIIR] Error: _eigenValues[i] < 0");
      return;
    }

    // NOTE: In our implementation, U^TMU = I, which is the same with Changxi's course notes, but different from the DyRT paper and WaveBlender's implementation
    // double omega = std::sqrt(_eigenValues[i] / _material.density);
    double omega = std::sqrt(_eigenValues[i] / 1.0 );
    double omega_di = _material.omega_di(omega);

    double xi = _material.xi(omega);
    double _epsilon = std::exp(-xi * omega * _dt);
    double _theta = omega_di * _dt;
    double _gamma = std::asin(xi);

    _coeff_qNew[i] = 2. * _epsilon * std::cos(_theta);
    _coeff_qOld[i] = _epsilon * _epsilon;
    _coeff_Q[i] = 2. / (3. * omega * omega_di) *
                  (_epsilon * std::cos(_theta + _gamma) -
                   _epsilon * _epsilon * std::cos(2. * _theta + _gamma)) /
                  _material.density;

    // std::cout << "coeff_qNew[" << i << "]: " << _coeff_qNew[i] << std::endl;
  }
}

// formula reference: https://www.cs.cornell.edu/~djames/papers/DyRT.pdf
double ModalIntegratorIIR::step(
    double time,
    const std::vector<Eigen::VectorXd> &forces,
    const std::vector<std::pair<int, int>> &objectIdx_and_vertexIdx)
{

  Eigen::VectorXd modalForce = Eigen::VectorXd::Zero(_N_modes); // force in modal space

  // Compute force at time step
  for (int i = 0; i < forces.size(); i++)
  {
    Eigen::VectorXd force = forces[i];
    auto [objIdx, vertexIdx] = objectIdx_and_vertexIdx[i];

    auto subforce = _eigenVectors.row(3 * vertexIdx).transpose() * force[0] +
                    _eigenVectors.row(3 * vertexIdx + 1).transpose() * force[1] +
                    _eigenVectors.row(3 * vertexIdx + 2).transpose() * force[2];

    modalForce += subforce.transpose();
  }

  // Step the system
  _q_next_next = _coeff_qNew.cwiseProduct(_q_next) -
                 _coeff_qOld.cwiseProduct(_q_curr) +
                 _coeff_Q.cwiseProduct(modalForce);

  _qDot_curr_plus = (_q_next - _q_curr) / _dt;                     // velocity
  _qDDot_curr = (_q_next + _q_prev - 2.0 * _q_curr) / (_dt * _dt); // acceleration

  // Update pointers
  _q_prev = _q_curr;
  _q_curr = _q_next;
  _q_next = _q_next_next;

  // std::cout << "modalForce: " << modalForce.transpose() << std::endl;
  // std::cout << "q_curr: " << _q_curr.transpose() << std::endl;

  // return _qDDot_curr.sum();
  return _q_next_next.sum();
}
