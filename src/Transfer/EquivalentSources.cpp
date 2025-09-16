#include <Transfer/SphericalHarmonics.h>
#include <Transfer/EquivalentSources.h>

#include <Eigen/Core>
#include <Eigen/SVD>

#include <iostream>
#include <complex>
#include <cmath>
#include "utils.h"

// #include "/opt/homebrew/Cellar/libomp/20.1.8/include/omp.h"

using namespace std::complex_literals; 
using namespace Transfer;

std::vector<Eigen::VectorXcd> EquivalentSources::computeTransferCoefficients(
      const Eigen::MatrixXd &V,
      const Eigen::MatrixXd &N,
      const Eigen::MatrixXd &sources,
      const std::vector<Eigen::MatrixXd> &displacements,
      Eigen::VectorXd omega_squareds)
{
    Utils::printInfo("[2] Solving for Equivalent Source coefficients");

    int num_modes          = omega_squareds.size();
    int num_sources        = sources.rows();
    double densityOfMedium = air_density;

    std::cout << "Number of sources: " << num_sources << std::endl;

    std::vector<Eigen::VectorXcd> cs;
    cs.reserve(num_modes);

    // get object diameter 
    double L = (V.colwise().maxCoeff() - V.colwise().minCoeff()).maxCoeff();

    // precompute angular part of spherial harmonics
    TIC(PrecomputeAngularCaches);
    std::vector<std::array<AngularCache,4>> sourceCaches;
    sourceCaches.reserve(num_sources);
    for (int s = 0; s < num_sources; ++s) {
        Eigen::Vector3d source       = sources.row(s).transpose();
        Eigen::MatrixXd SphWRTSource = SphericalHarmonics::GetSphericalCoordinates(V, source);
        Eigen::MatrixXd N_Sph        = SphericalHarmonics::ConvertNormalToLocalSpherical(V, N, source);

        std::array<AngularCache,4> caches = {
            SphericalHarmonics::PrecomputeAngular(0,  0, SphWRTSource, N_Sph),
            SphericalHarmonics::PrecomputeAngular(1, -1, SphWRTSource, N_Sph),
            SphericalHarmonics::PrecomputeAngular(1,  0, SphWRTSource, N_Sph),
            SphericalHarmonics::PrecomputeAngular(1,  1, SphWRTSource, N_Sph)
        };
        sourceCaches.push_back(std::move(caches));
    }
    TOC(PrecomputeAngularCaches);

    // // method of known solutions
    // Eigen::Vector3d sample_source;
    //     sample_source << 0, 0, 0;
    // Eigen::MatrixXd Sph_V_wrt_0 = SphericalHarmonics::GetSphericalCoordinates(V, sample_source);
    // Eigen::MatrixXd N_Origin_Spherical = SphericalHarmonics::ConvertNormalToLocalSpherical(V, N, Eigen::Vector3d::Zero());
    // AngularCache cache_monopole = SphericalHarmonics::PrecomputeAngular(1,  0, Sph_V_wrt_0, N_Origin_Spherical);

    for (int mode = 0; mode < num_modes; ++mode)
    {
        std::cout << "--------------------" << std::endl;

        double omega                = sqrt(omega_squareds[mode]);
        double k                    = omega / speed_of_sound;
        const Eigen::MatrixXd &disp = displacements[mode];

        // compute normal displacements
        Eigen::VectorXd u_n = (disp.array() * N.array()).rowwise().sum();

        // build the derivative matrix 
        TIC(BuildDerivativeMatrix);
        Eigen::MatrixXcd A(V.rows(), 4 * num_sources);
        for (int s = 0; s < num_sources; ++s)
        {
            int offset = 4 * s;
            const auto &caches = sourceCaches[s];

            A.col(offset + 0) = SphericalHarmonics::EvaluateRadial(k, caches[0]); // (0, 0)
            A.col(offset + 1) = SphericalHarmonics::EvaluateRadial(k, caches[1]); // (1,-1)
            A.col(offset + 2) = SphericalHarmonics::EvaluateRadial(k, caches[2]); // (1, 0)
            A.col(offset + 3) = SphericalHarmonics::EvaluateRadial(k, caches[3]); // (1, 1)
        }
        TOC(BuildDerivativeMatrix);

        // define the right hand side boundary condition of linear system 
        Eigen::VectorXcd b = densityOfMedium * omega * omega * u_n.cast<std::complex<double>>();

        // // analytic monopole BC for testing
        // Eigen::VectorXcd b = SphericalHarmonics::EvaluateRadial(k, cache_monopole);

        TIC(SolveLinearSystem);

        // this seems fine
        double epsilon       = 1e-2 * A.norm();
        Eigen::MatrixXcd AtA = A.adjoint() * A + epsilon * epsilon * Eigen::MatrixXcd::Identity(A.cols(), A.cols());
        Eigen::VectorXcd Atb = A.adjoint() * b;
        Eigen::VectorXcd c   = AtA.ldlt().solve(Atb);

        TOC(SolveLinearSystem);


        Eigen::VectorXcd residual = A * c - b;

        // print stuff
        std::cout << "Mode: " << mode << std::endl;
        std::cout << "Linear Frequency: " << omega / (2 * M_PI) << " Hz"<< std::endl;
        std::cout << "kL: " << k*L << std::endl;
        std::cout << "c norm: " << c.norm() << std::endl;
        std::cout << "Relative residual: " << residual.norm() / b.norm() << std::endl;

        cs.push_back(std::move(c));
    }
    
    return cs;
}

