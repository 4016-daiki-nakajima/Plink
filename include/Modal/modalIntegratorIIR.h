#pragma once

#include "modalIntegratorBase.h"

namespace SimpleModal
{

  class ModalIntegratorIIR : public ModalIntegratorBase
  {
  public:
    ModalIntegratorIIR(const double dt, const Eigen::VectorXd &eigenValues,
                       const Eigen::MatrixXd &eigenVectors,
                       const ModalMaterial &material);

    double step(double time,
                const std::vector<Eigen::VectorXd> &forces,
                const std::vector<std::pair<int, int>> &objectIdx_and_vertexIdx) override;

    inline Eigen::VectorXd get_q_current() const { return _q_curr; }

  protected:
    int _N_modes;
    Eigen::VectorXd _q_prev;
    Eigen::VectorXd _q_curr;
    Eigen::VectorXd _q_next;
    Eigen::VectorXd _q_next_next;

    Eigen::VectorXd _qDot_curr_plus; // modal velocity
    Eigen::VectorXd _qDDot_curr;     // modal acceleration

    Eigen::VectorXd _coeff_qNew;
    Eigen::VectorXd _coeff_qOld;
    Eigen::VectorXd _coeff_Q;
  };
}; // namespace SimpleModal
