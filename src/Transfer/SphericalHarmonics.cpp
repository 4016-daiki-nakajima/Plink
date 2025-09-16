#include <Transfer/SphericalHarmonics.h>

#include <Eigen/Core>

#include <iostream>
#include <complex>
#include <cmath>

// to use these must uncomment the "find_package(Boost REQUIRED)"" in CMakeLists.txt
// #include <boost/math/special_functions/bessel.hpp>
// #include <boost/math/special_functions/legendre.hpp>
// #include <boost/math/special_functions/spherical_harmonic.hpp>

using namespace std::complex_literals; 
using namespace Transfer;

Eigen::Vector3d SphericalHarmonics::get_spherical_coords(const Eigen::Vector3d &cartesian)
{
    double r     = cartesian.norm();
    double theta = std::acos(cartesian.y() / r);
    double phi   = std::atan2(cartesian.z(), cartesian.x());
    
    return Eigen::Vector3d(r, theta, phi);
}

Eigen::Vector3d SphericalHarmonics::get_cartesian_coords(const Eigen::Vector3d &spherical)
{
    double r     = spherical(0);
    double theta = spherical(1);
    double phi   = spherical(2);

    double x = r * sin(theta) * cos(phi);
    double y = r * cos(theta);
    double z = r * sin(theta) * sin(phi);

    return Eigen::Vector3d(x, y, z);
}


Eigen::MatrixXd SphericalHarmonics::GetSphericalCoordinates(const Eigen::MatrixXd& V, const Eigen::Vector3d &source) {

    Eigen::MatrixXd spherical_coordinates(V.rows(), 3);
    for (int i = 0; i < V.rows(); ++i) {
        Eigen::Vector3d shifted = V.row(i).transpose() - source;
        spherical_coordinates.row(i) = Transfer::SphericalHarmonics::get_spherical_coords(shifted).transpose();
    }
    return spherical_coordinates;
}

static Eigen::MatrixXd GetCartesianCoordinates(const Eigen::MatrixXd& Spherical) {

    Eigen::MatrixXd cartesian_coordinates(Spherical.rows(), 3);
    for (int i = 0; i < Spherical.rows(); ++i) {
        cartesian_coordinates.row(i) = Transfer::SphericalHarmonics::get_cartesian_coords(Spherical.row(i).transpose()).transpose();
    }
    return cartesian_coordinates;
}

static Eigen::Vector3d convert_normal_to_spherical(const Eigen::Vector3d &origin, const Eigen::Vector3d &N)
{
    Eigen::Vector3d spherical_coords = Transfer::SphericalHarmonics::get_spherical_coords(origin);
    if (abs(spherical_coords(0)) < 1e-6) return Eigen::Vector3d(0, 0, 0); // handle singularity

    double r     = spherical_coords(0);
    double theta = spherical_coords(1);
    double phi  = spherical_coords(2);

    double sin_theta = sin(theta);
    double cos_theta = cos(theta);
    double sin_phi   = sin(phi);
    double cos_phi   = cos(phi);

    // Basis vectors for y-up spherical coordinates
    Eigen::Vector3d e_r    ( sin_theta*cos_phi,  cos_theta,  sin_theta*sin_phi );
    Eigen::Vector3d e_theta( cos_theta*cos_phi, -sin_theta,  cos_theta*sin_phi );
    Eigen::Vector3d e_phi  ( -sin_phi,           0,          cos_phi           );

    // normalize N
    Eigen::Vector3d N_normalized = N.normalized();

    // project N onto spherical coordinates
    Eigen::Vector3d N_spherical;
    N_spherical(0) = N_normalized.dot(e_r);
    N_spherical(1) = N_normalized.dot(e_theta);
    N_spherical(2) = N_normalized.dot(e_phi);

    return N_spherical;
}

Eigen::MatrixXd SphericalHarmonics::ConvertNormalToLocalSpherical(
    const Eigen::MatrixXd& V, 
    const Eigen::MatrixXd& N, 
    const Eigen::Vector3d &source) 
{
    
    Eigen::MatrixXd spherical(V.rows(), 3);
    for (int i = 0; i < V.rows(); ++i) {
        Eigen::Vector3d rel = V.row(i).transpose() - source;
        spherical.row(i) = convert_normal_to_spherical(rel, N.row(i)).transpose();
    }
    return spherical;
}

inline double factorial(int x) {
    if (x < 0) {
        throw std::domain_error("Factorial not defined for negative numbers");
    }
    return std::tgamma(x + 1);
}

inline double NormalizationFactor(int ell, int m) {
    return sqrt(((2.0 * ell + 1.0) / (4.0 * M_PI)) * ( factorial(ell - m) / factorial(ell + m) ));
}

static Eigen::ArrayXcd hankel_2(int ell, const Eigen::ArrayXd &kr) {
    Eigen::ArrayXcd result(kr.size());
    if (ell == 0) { // h0^(2) = ie^{-i kr}/kr
        result = 1i * (-1i * kr).exp() / kr.cast<std::complex<double>>();  
    }
    else if (ell == 1) {  // h1^(2) = - e^{-i kr} * (kr - i)/kr**2
        result = -(kr.cast<std::complex<double>>() - 1i) * (-1i * kr).exp() / (kr * kr).cast<std::complex<double>>(); 
    }
    else {
        std::cerr << "Not implemented for ell: " << ell << std::endl;
        result.setZero();
    }
    return result;
}


static Eigen::ArrayXcd DV_hankel_2(int ell, const Eigen::ArrayXd &kr) {
    Eigen::ArrayXcd result(kr.size());

    if (ell == 0) {
        result = (-(1i * kr)).exp() * (kr - 1i) / (kr * kr).cast<std::complex<double>>(); // h0^(2) = ie^{-i kr}/kr
    } else if (ell == 1) {
        result = (-(1i * kr)).exp() *
                 ((1i * kr * kr) + (2 * kr) - 2i) / (kr * kr * kr).cast<std::complex<double>>();  // h1^(2) = - e^{-i kr} * (kr - i)/kr**2
    } else {
        std::cerr << "Not implemented for ell: " << ell << std::endl;
        result.setZero();
    }

    return result;
}

static Eigen::ArrayXd P_nm(int ell, int m, const Eigen::ArrayXd &theta) {
    Eigen::ArrayXd result(theta.size());
    if (ell == 0) {
        result.setOnes();
    } else if (ell == 1) {
        if      (m == -1) result = 0.5 * theta.sin();
        else if (m ==  0) result = theta.cos();
        else if (m ==  1) result = -theta.sin();
    } else {
        std::cerr << "Not implemented for ell: " << ell << ", m: " << m << std::endl;
        result.setZero();
    }
    return result;
}

static Eigen::ArrayXd DV_P_nm(int ell, int m, const Eigen::ArrayXd &theta) {
    Eigen::ArrayXd result(theta.size());
    if (ell == 0) {
        result.setZero();
    } else if (ell == 1) {
        if      (m == -1) result = 0.5 * theta.cos();
        else if (m ==  0) result = -theta.sin();
        else if (m ==  1) result = -theta.cos();
    } else {
        std::cerr << "Not a dipole! ell: " << ell << ", m: " << m << std::endl;
        result.setZero();
    }
    return result;
}

AngularCache SphericalHarmonics::PrecomputeAngular(
    int ell, int m,
    const Eigen::MatrixXd &SphericalCoordsOfVwithRespectToSource, 
    const Eigen::MatrixXd &NormalInLocalSphericalCoordinates)     
{
    AngularCache cache;
    cache.ell = ell;
    cache.m   = m;
    cache.NormalizationFactor = NormalizationFactor(ell, m);

    cache.r              = SphericalCoordsOfVwithRespectToSource.col(0).array();
    Eigen::ArrayXd theta = SphericalCoordsOfVwithRespectToSource.col(1).array();
    Eigen::ArrayXd phi   = SphericalCoordsOfVwithRespectToSource.col(2).array();

    cache.sin_theta = theta.sin();
    cache.exp_imphi = (1i * static_cast<double>(m) * phi).exp();

    cache.Pnm  = P_nm(ell, m, theta);
    cache.dPnm = DV_P_nm(ell, m, theta);

    cache.N_spherical = NormalInLocalSphericalCoordinates;

    return cache;
}


Eigen::VectorXcd SphericalHarmonics::EvaluateRadial(
    double k,
    const AngularCache &cache)
{
    Eigen::ArrayXd  kr         = k * cache.r;

    Eigen::ArrayXcd hankel     = hankel_2(cache.ell, kr);
    Eigen::ArrayXcd dhankel    = DV_hankel_2(cache.ell, kr);

    Eigen::ArrayXcd grad_r     = k * dhankel * cache.Pnm;
    Eigen::ArrayXcd grad_theta = (1.0 / cache.r) * hankel * cache.dPnm;
    Eigen::ArrayXcd grad_phi   = ((std::complex<double>(0, cache.m) / (cache.r * cache.sin_theta)) 
                                   * hankel * cache.Pnm );

    grad_phi = (cache.sin_theta.abs() < 1e-6).select(std::complex<double>(0,0), grad_phi);

    int N = cache.r.size();
    Eigen::VectorXcd results(N);

    for (int i = 0; i < N; ++i) {
        Eigen::Vector3cd g;
        g << grad_r(i), grad_theta(i), grad_phi(i);

        Eigen::Vector3d N_s = cache.N_spherical.row(i);
        results(i) = g.dot(N_s.cast<std::complex<double>>()) * cache.exp_imphi(i);
    }

    return cache.NormalizationFactor * results;
}

std::complex<double> SphericalHarmonics::Psi(
    int ell, int m, double k,
    const Eigen::Vector3d &x,
    const Eigen::Vector3d &x0)
{
    Eigen::Vector3d spherical_coords = Transfer::SphericalHarmonics::get_spherical_coords(x - x0);
    double r     = spherical_coords(0);
    double theta = spherical_coords(1);
    double phi   = spherical_coords(2);
    double kr    = k * r;

    // wrap scalars into length-1 arrays since the functions are defined for arrays
    Eigen::ArrayXd thetaArr(1);
    Eigen::ArrayXd krArr(1);
    Eigen::ArrayXd phiArr(1);

    thetaArr << theta;
    krArr    << kr;
    phiArr   << phi;

    Eigen::ArrayXcd result = hankel_2(ell, krArr) * NormalizationFactor(ell, m) * P_nm(ell, m, thetaArr) * (1i * static_cast<double>(m) * phiArr).exp();
    
    return result(0);
}


Eigen::ArrayXcd SphericalHarmonics::Psi(
    int ell, int m, double k,
    const Eigen::MatrixXd &SphericalCoordinatesWRTSource)
{
    Eigen::Map<const Eigen::ArrayXd> r(SphericalCoordinatesWRTSource.col(0).data(), SphericalCoordinatesWRTSource.rows());
    Eigen::Map<const Eigen::ArrayXd> theta(SphericalCoordinatesWRTSource.col(1).data(), SphericalCoordinatesWRTSource.rows());
    Eigen::Map<const Eigen::ArrayXd> phi(SphericalCoordinatesWRTSource.col(2).data(), SphericalCoordinatesWRTSource.rows());

    Eigen::ArrayXcd PsiArr =   hankel_2(ell, k * r) 
                             * NormalizationFactor(ell, m) 
                             * P_nm(ell, m, theta) 
                             * (1i * static_cast<double>(m) * phi).exp();

    return PsiArr;
}


Eigen::ArrayXcd SphericalHarmonics::Pressure(
    double k,
    const Eigen::VectorXcd &coeffs,
    const Eigen::MatrixXd &V,       
    const Eigen::MatrixXd &sources) 
{
    int num_points  = V.rows();
    int num_sources = sources.rows();

    Eigen::ArrayXcd pressure = Eigen::ArrayXcd::Zero(num_points);

    for (int s = 0; s < num_sources; ++s) {
        int offset = 4 * s;
        Eigen::MatrixXd sph = SphericalHarmonics::GetSphericalCoordinates(V, sources.row(s));

        pressure += SphericalHarmonics::Psi(0,  0, k, sph) * coeffs[offset + 0];
        pressure += SphericalHarmonics::Psi(1, -1, k, sph) * coeffs[offset + 1];
        pressure += SphericalHarmonics::Psi(1,  0, k, sph) * coeffs[offset + 2];
        pressure += SphericalHarmonics::Psi(1,  1, k, sph) * coeffs[offset + 3];
    }

    return pressure;
}







// std::complex<double> SphericalHarmonics::Pressure(
//     double k,
//     const Eigen::Vector3d  &x,
//     const Eigen::VectorXcd &coeffs,
//     const Eigen::MatrixXd  &sources)
// {
//     std::complex<double> pressure = 0.0;

//     int num_sources = sources.rows();

//     for (int s = 0; s < num_sources; ++s)
//     {
//         Eigen::Vector3d source = sources.row(s);
//         int offset = 4 * s;

//         pressure += SphericalHarmonics::Psi(0,  0, k, x, source) * coeffs[offset + 0];
//         pressure += SphericalHarmonics::Psi(1, -1, k, x, source) * coeffs[offset + 1];
//         pressure += SphericalHarmonics::Psi(1,  0, k, x, source) * coeffs[offset + 2];
//         pressure += SphericalHarmonics::Psi(1,  1, k, x, source) * coeffs[offset + 3];
//     }

//     return pressure;
// }






// // try computing gradient in spherical coordinates, convert BACK to cartesian, and then dot product
// // this way we can keep everything in global coordinates
// Eigen::VectorXcd SphericalHarmonics::NormalDVPsiGradient(
//     int ell, int m, double k,
//     const Eigen::MatrixXd &V, 
//     const Eigen::MatrixXd &N,
//     const Eigen::Vector3d &source)
// {
//     // Spherical coords of (x - source)
//     Eigen::MatrixXd Sph = GetSphericalCoordinates(V, source);

//     Eigen::ArrayXd r     = Sph.col(0).array();
//     Eigen::ArrayXd theta = Sph.col(1).array();
//     Eigen::ArrayXd phi   = Sph.col(2).array();

//     // Guard tiny r to avoid 1/r blow-ups
//     const double r_eps = 1e-12;
//     r = r.max(r_eps);

//     Eigen::ArrayXd  sin_theta = theta.sin();
//     Eigen::ArrayXcd exp_imphi = (1i * static_cast<double>(m) * phi).exp();

//     Eigen::ArrayXd kr = k * r;

//     Eigen::ArrayXcd Pnm  = P_nm(ell, m, theta);
//     Eigen::ArrayXcd dPnm = DV_P_nm(ell, m, theta);
//     Eigen::ArrayXcd h2   = hankel_2(ell, kr);
//     Eigen::ArrayXcd dh2  = DV_hankel_2(ell, kr);

//     // Spherical gradient components (basis is {e_r, e_theta, e_phi} at x - source)
//     Eigen::ArrayXcd grad_r     = k * dh2 * Pnm;
//     Eigen::ArrayXcd grad_theta = (1.0 / r) * h2 * dPnm;
//     Eigen::ArrayXcd grad_phi   = ( (std::complex<double>(0.0, static_cast<double>(m)) / (r * sin_theta)) * h2 * Pnm );

//     // Handle polar axis where sin(theta) ~ 0
//     grad_phi = (sin_theta.abs() < 1e-9).select(std::complex<double>(0.0, 0.0), grad_phi);

//     const int n = V.rows();
//     Eigen::VectorXcd out(n);

//     for (int i = 0; i < n; ++i) {
//         const double th = theta(i);
//         const double ph = phi(i);

//         const double st = std::sin(th), ct = std::cos(th);
//         const double sp = std::sin(ph), cp = std::cos(ph);

//         // y-up spherical basis at the field point (using your conventions)
//         const Eigen::Vector3d e_r    ( st*cp,  ct,  st*sp );
//         const Eigen::Vector3d e_theta( ct*cp, -st,  ct*sp );
//         const Eigen::Vector3d e_phi  ( -sp,    0,   cp   );

//         // Rotate spherical components to Cartesian
//         const Eigen::Vector3cd grad_cart =
//               grad_r(i)     * e_r
//             + grad_theta(i) * e_theta
//             + grad_phi(i)   * e_phi;

//         // Use the **global** unit normal at the point — do NOT subtract the source
//         Eigen::Vector3d n_hat = N.row(i);
//         const double n_norm = n_hat.norm();
//         if (n_norm > 0.0) n_hat /= n_norm; else n_hat.setZero();

//         out(i) = grad_cart.dot(n_hat.cast<std::complex<double>>()) * exp_imphi(i);
//     }
//     return NormalizationFactor(ell, m) * out;
// }





