
#include <Eigen/Core>
#include <Eigen/Eigenvalues>

#include <iostream>
#include <cmath>
#include "Utils/utils.h"

#include <Transfer/Power.h>
#include <Transfer/SphericalHarmonics.h>
// #include <Transfer/EquivalentSources.h>


using namespace std::complex_literals; 
using namespace Transfer;

Eigen::VectorXd Power::computePowers(
      const std::vector<Eigen::VectorXcd>& coefficients,
      const Eigen::VectorXd& omega_squareds)
{
    int num_modes = coefficients.size();
    Eigen::VectorXd powers(num_modes);

    double factor = 1 / (speed_of_sound * air_density);
    Eigen::VectorXd k_squareds = omega_squareds / (speed_of_sound * speed_of_sound);

    Eigen::VectorXd coefficient_norms = Eigen::VectorXd::Zero(num_modes);
    for (int mode = 0; mode < num_modes; ++mode) {
        const Eigen::VectorXcd& cs = coefficients[mode];
        int num_sources = cs.size() / 4;

        Eigen::Map<const Eigen::MatrixXcd> coeff_matrix(cs.data(), 4, num_sources);
        Eigen::VectorXcd source_sums = coeff_matrix.colwise().sum();
        double norm_sum = source_sums.squaredNorm();

        powers[mode] = factor * k_squareds[mode] * norm_sum;
    }

    return powers;
}





static Eigen::MatrixXd fibonacci_sphere(double R, int num_points)
{
    Eigen::MatrixXd points(num_points, 3);
    double golden_angle = M_PI * (3. - sqrt(5.)); // golden angle in radians
    for (int i = 0; i < num_points; ++i) {
        double y = 1 - (i / (double)(num_points - 1)) * 2; // y goes from 1 to -1
        double r = sqrt(1 - y * y); // radius at y

        double theta = golden_angle * i; 
        double x     = cos(theta) * r;
        double z     = sin(theta) * r;

        points(i, 0) = x;
        points(i, 1) = y;
        points(i, 2) = z;
    }
    return R * points;
}

double capIntegral(const Eigen::Vector3d& r,
                   const Eigen::MatrixXd& samples,
                   const Eigen::VectorXd& s_vals,
                   double cap_angle_rad)
{
    double cos_th = std::cos(cap_angle_rad);
    double sum    = 0.0;
    for (int i = 0; i < samples.rows(); ++i) {
        Eigen::Vector3d n = samples.row(i);
        if (std::abs(r.dot(n)) >= cos_th) {
            sum += s_vals[i];
        }
    }
    return sum;
}


std::vector<Eigen::Vector3d> Power::computeEquipowerDipoleDirections(
    const Eigen::MatrixXd &sources,
    const std::vector<Eigen::VectorXcd>& coefficients,
    const Eigen::VectorXd& omega_squareds)
{
    Utils::printInfo("[3] Equipower Dipole Directions");

    int num_modes = coefficients.size();
    std::vector<Eigen::Vector3d> dipole_directions(num_modes);

    // use fibonacci sphere to uniformly sample points on a sphere
    double R = 10.0;
    int num_samples = 1000;
    Eigen::MatrixXd samples_on_sphere = fibonacci_sphere(R, num_samples);

    // Precompute sin(theta) weights for spherical integration
    Eigen::ArrayXd sin_theta(num_samples);
    for (int i = 0; i < num_samples; ++i) {
        Eigen::Vector3d n = samples_on_sphere.row(i);
        double theta      = acos(n.normalized().z()); 
        sin_theta[i]      = sin(theta);
    }

    for (int mode = 0; mode < num_modes; ++mode) {
        const Eigen::VectorXcd &c = coefficients[mode];
        double k = sqrt(omega_squareds[mode]) / speed_of_sound;

        // compute pressures on the sphere and define the symmetric function S
        Eigen::ArrayXcd P1 = SphericalHarmonics::Pressure(k, c,  samples_on_sphere, sources);      
        Eigen::ArrayXcd P2 = SphericalHarmonics::Pressure(k, c, -samples_on_sphere, sources);  
        Eigen::ArrayXd  S  = P1.cwiseAbs2() + P2.cwiseAbs2();  

        // integrate to get the dipole moment
        Eigen::ArrayXd  w = S * sin_theta; // integration weights
        Eigen::Matrix3d D = samples_on_sphere.transpose() * w.matrix().asDiagonal() * samples_on_sphere;

        // solve eigenproblem
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eig(D);
        Eigen::VectorXd evals = eig.eigenvalues();
        Eigen::Matrix3d evecs = eig.eigenvectors();

        // index of the largest eigenvalue
        int idx_max = 0;
        evals.normalize();
        evals.maxCoeff(&idx_max);

        // check if any eigenvalues are nearly equal
        double tol = 1e-2;
        bool nearly_equal = ( (std::abs(evals(0) - evals(1)) < tol) ||
                              (std::abs(evals(0) - evals(2)) < tol) ||
                              (std::abs(evals(1) - evals(2)) < tol) );

        // if yes, then the principal eigenvector is ambiguous
        if (nearly_equal) {
            // do cap integral check for each eigenvector
            double cap_angle = M_PI / 6.0; // small angle around the eigenvector
            std::vector<double> cap_vals(3);
            for (int j = 0; j < 3; ++j) {
                Eigen::Vector3d u = evecs.col(j).normalized();
                cap_vals[j]       = capIntegral(u, samples_on_sphere, S, cap_angle);
            }
            // take index of max cap value
            idx_max = std::distance(cap_vals.begin(), std::max_element(cap_vals.begin(), cap_vals.end()));
        }
        dipole_directions[mode] = evecs.col(idx_max).normalized();
    }

    return dipole_directions;
}
