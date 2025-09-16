#pragma once

#include <cmath>
#include <Eigen/Core>


struct AngularCache {
    int ell, m;
    double NormalizationFactor;
    Eigen::ArrayXcd Pnm;
    Eigen::ArrayXcd dPnm;
    Eigen::ArrayXcd exp_imphi;
    Eigen::ArrayXd  r;
    Eigen::ArrayXd  sin_theta;
    Eigen::MatrixXd N_spherical;
};


namespace Transfer
{
  class SphericalHarmonics
  {
  public:

    // ---------------------------------------------------------------------
    // coordinate conversions
    // ---------------------------------------------------------------------

    static Eigen::Vector3d get_spherical_coords(
        const Eigen::Vector3d &cartesian);

    static Eigen::Vector3d get_cartesian_coords(
        const Eigen::Vector3d &spherical);

    static Eigen::MatrixXd GetSphericalCoordinates(
        const Eigen::MatrixXd& V, 
        const Eigen::Vector3d &source);

    static Eigen::MatrixXd ConvertNormalToLocalSpherical(
        const Eigen::MatrixXd& V, 
        const Eigen::MatrixXd& N, 
        const Eigen::Vector3d &source);

    // ---------------------------------------------------------------------
    // Sphereical Harmonic Bases
    // ---------------------------------------------------------------------

    static std::complex<double> Psi(
        int ell, int m, double k,
        const Eigen::Vector3d &x,
        const Eigen::Vector3d &x0);

    static Eigen::ArrayXcd Psi(
        int ell, int m, double k,
        const Eigen::MatrixXd &SphericalCoordinatesWRTSource);

    static Eigen::VectorXcd NormalDVPsiGradient(
        int ell, int m, double k,
        const Eigen::MatrixXd &V, 
        const Eigen::MatrixXd &N,
        const Eigen::Vector3d &source);

    static AngularCache PrecomputeAngular(
        int ell, int m,
        const Eigen::MatrixXd &SphericalCoordsOfVwithRespectToSource,
        const Eigen::MatrixXd &NormalInLocalSphericalCoordinates);

    static Eigen::VectorXcd EvaluateRadial(
        double k,
        const AngularCache &cache);



    // ---------------------------------------------------------------------
    // Pressure
    // ---------------------------------------------------------------------

    static Eigen::ArrayXcd Pressure(
        double k,
        const Eigen::VectorXcd &coeffs,
        const Eigen::MatrixXd  &V,
        const Eigen::MatrixXd  &sources);
  };

    
} // namespace Transfer