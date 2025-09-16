
#include <Eigen/Core>
#include <Eigen/LU>

#include <iostream>
#include "utils.h"

#include <Geometry/BEM.h>

using namespace Plink;

Eigen::Matrix3d BEM::ComputeAddedMass(
  const Eigen::MatrixXd &V,
  const Eigen::MatrixXi &F,
  const Eigen::MatrixXd &per_face_normals,
  const Eigen::VectorXd &per_face_areas)
{
    // exterior direct BEM
    int num_faces = F.rows();

    // exterior BVP requires normals pointing inwards
    Eigen::MatrixXd flipped_normals = -per_face_normals;

    // compute centroids of faces
    Eigen::MatrixXd centroids(num_faces, 3);
    for (int i = 0; i < num_faces; ++i) {
        centroids.row(i) = (V.row(F(i, 0)) + V.row(F(i, 1)) + V.row(F(i, 2))) / 3.0;
    }

    // preallocate H and G matrices
    Eigen::MatrixXd H = Eigen::MatrixXd::Zero(num_faces, num_faces);
    Eigen::MatrixXd G = Eigen::MatrixXd::Zero(num_faces, num_faces);

    // preallocate B matrix for Neumann BC
    Eigen::MatrixXd B = Eigen::MatrixXd::Zero(num_faces, 3);

    // build H and G matrices
    for (int i = 0; i < num_faces; ++i) {
        Eigen::Vector3d xi = centroids.row(i);
        for (int j = 0; j < num_faces; ++j) {
            Eigen::Vector3d yj = centroids.row(j);
            Eigen::Vector3d nj = flipped_normals.row(j).normalized();
            double aj = per_face_areas(j);

            double Gval;
            double Hval;
            if (i == j) { // diagonal terms
                Gval = 2.0 * std::sqrt(aj) / (3.0 * 4.0 * M_PI);
                Hval = -0.5;
            } else {
                Eigen::Vector3d r = xi - yj;
                double rlen = r.norm();

                Gval =  aj / (4.0 * M_PI * rlen);
                Hval = -aj * (r.dot(nj)) / (4.0 * M_PI * rlen * rlen * rlen);
            }
            G(i,j) = Gval;
            H(i,j) = Hval;

            // accumulate b_i
            B.row(i) -= (Gval * nj).transpose();
        }
    }

    // pre-factor H
    Eigen::PartialPivLU<Eigen::MatrixXd> Hsolver;
    Hsolver.compute(H);

    Eigen::MatrixXd U(num_faces, 3);
    for (int k = 0; k < 3; ++k) {
        Eigen::VectorXd rhs = B.col(k);
        U.col(k) = Hsolver.solve(rhs);
    }

    // accumulate added mass
    Eigen::Matrix3d M = Eigen::Matrix3d::Zero();
    for (int i = 0; i < num_faces; ++i) {
        double a = per_face_areas(i);
        Eigen::Vector3d n = flipped_normals.row(i).normalized(); 
        Eigen::Vector3d u = U.row(i).transpose();                  
        M += u * n.transpose() * a;                                
    }

    return M;
}
