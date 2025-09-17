#include "Geometry/FEM.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "Utils/utils.h"

using namespace SimpleModal;
// Helper function to compute the volume of a tetrahedron
double FEM::computeTetrahedronVolume(const Eigen::MatrixXd &TV, const Eigen::Vector4i &tet)
{
    Eigen::Matrix<double, 3, 3> Dm;
    Dm.col(0) = TV.row(tet[1]) - TV.row(tet[0]);
    Dm.col(1) = TV.row(tet[2]) - TV.row(tet[0]);
    Dm.col(2) = TV.row(tet[3]) - TV.row(tet[0]);
    return std::abs(Dm.determinant()) / 6.0;
}


Eigen::Matrix<double, 6, 6> computeElasticityMatrix(double E, double nu)
{
    // Construct the material stiffness matrix
    Eigen::Matrix<double, 6, 6> D;
    D.setZero();
    D.block<3, 3>(0, 0) << 1 - nu, nu, nu,
        nu, 1 - nu, nu,
        nu, nu, 1 - nu;
    D.block<3, 3>(3, 3) = Eigen::Matrix3d::Identity() * (0.5 - nu);
    D *= E / (1 + nu) / (1 - 2 * nu);
    // std::cout << Utils::CYAN << "D: " << D << Utils::RESET << std::endl;

    return D;
}

// Helper function to compute the stiffness matrix for a single tetrahedron
Eigen::Matrix<double, 12, 12> computeElementStiffnessMatrix(const Eigen::MatrixXd &TV, const Eigen::Vector4i &tet, const Eigen::Matrix<double, 6, 6> &D)
{
    // Extract the vertices of the tetrahedron
    Eigen::Matrix<double, 3, 4> X;
    for (int i = 0; i < 4; ++i)
    {
        X.col(i) = TV.row(tet[i]);
    }

    // Compute the deformation gradient
    Eigen::Matrix<double, 3, 3> Dm;
    Dm.col(0) = X.col(1) - X.col(0);
    Dm.col(1) = X.col(2) - X.col(0);
    Dm.col(2) = X.col(3) - X.col(0);
    // std::cout << Utils::YELLOW << "Dm: " << Dm << Utils::RESET << std::endl;

    double volume = std::abs(Dm.determinant()) / 6.0;
    Eigen::Matrix<double, 3, 3> Dm_inv = Dm.inverse();
    // std::cout << Utils::CYAN << "Dm_inv: " << Dm_inv << Utils::RESET << std::endl;

    // Compute the B matrix
    Eigen::Matrix<double, 6, 12> B;
    B.setZero();

    auto delta = [](int i, int j) -> double
    { return i == j ? 1 : 0; };
    for (int i = 0; i < 4; ++i)
    {
        Eigen::Matrix<double, 3, 1> grad_N = Dm_inv.transpose() * Eigen::Matrix<double, 3, 1>(delta(i, 1) - delta(i, 0), delta(i, 2) - delta(i, 0), delta(i, 3) - delta(i, 0));
        //std::cout << Utils::YELLOW << "grad_N: " << grad_N << Utils::RESET << std::endl;
        B.block<6, 3>(0, 3 * i) <<  grad_N[0], 0, 0,
                                    0, grad_N[1], 0,
                                    0, 0, grad_N[2],
                                    grad_N[1], grad_N[0], 0,
                                    0, grad_N[2], grad_N[1],
                                    grad_N[2], 0, grad_N[0];
    }
    // std::cout << Utils::CYAN << "B: " << B << Utils::RESET << std::endl;

    // Construct the stiffness matrix for the tetrahedron
    Eigen::Matrix<double, 12, 12> K = Eigen::Matrix<double, 12, 12>::Zero();
    // Compute the element stiffness matrix
    K = volume * B.transpose() * D * B;

    return K;
}

Eigen::SparseMatrix<double> FEM::computeStiffnessMatrix(const TetMesh &mesh, const double &young_modulus, const double &poisson_ratio)
{
    const int numVertices = mesh.TV.rows();
    const int ndim = 3;
    Eigen::SparseMatrix<double> stiffnessMatrix(ndim * numVertices, ndim * numVertices);
    std::vector<Eigen::Triplet<double>> triplets;

    Eigen::Matrix<double, 6, 6> D = computeElasticityMatrix(young_modulus, poisson_ratio);

    for (int i = 0; i < mesh.TT.rows(); ++i) // i is the index of the tetrahedron
    {
        Eigen::Vector4i tet = mesh.TT.row(i);
        Eigen::Matrix<double, 12, 12> K = computeElementStiffnessMatrix(mesh.TV, tet, D);

        // Add contributions to the global stiffness matrix
        for (int j = 0; j < 4; ++j) // traverse the vertices of the tetrahedron
        {
            for (int k = 0; k < 4; ++k)
            {
                // traverse the 3x3 DoFs of one vertex
                for (int a = 0; a < ndim; ++a)
                {
                    for (int b = 0; b < ndim; ++b)
                    {
                        triplets.emplace_back(ndim * tet[j] + a, ndim * tet[k] + b, K(ndim * j + a, ndim * k + b));
                    }
                }
            }
        }
    }

    stiffnessMatrix.setFromTriplets(triplets.begin(), triplets.end());
    return stiffnessMatrix;
}

Eigen::SparseMatrix<double> FEM::computeMassMatrix(const TetMesh &mesh, const double &density)
{
    const int numVertices = mesh.TV.rows();
    const int ndim = 3;
    // 3N \times 3N lumped mass matrix
    Eigen::SparseMatrix<double> massMatrix(ndim * numVertices, ndim * numVertices);

    std::vector<Eigen::Triplet<double>> triplets;
    for (int i = 0; i < mesh.TT.rows(); ++i)
    {
        Eigen::Vector4i tet = mesh.TT.row(i);
        double volume = computeTetrahedronVolume(mesh.TV, tet);
        double mass = volume * density / 4.0;

        for (int j = 0; j < 4; ++j)
        {
            for (int k = 0; k < ndim; ++k) // [x, y, z]
            {
                triplets.emplace_back(ndim * tet[j] + k, ndim * tet[j] + k, mass);
            }
        }
    }

    massMatrix.setFromTriplets(triplets.begin(), triplets.end());
    return massMatrix;
}
