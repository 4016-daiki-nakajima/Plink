#include "Transfer/Sampling.h"
#include "FEM.h" 
#include <random>
#include <algorithm>
#include <igl/offset_surface.h>
#include <igl/random_points_on_mesh.h>


using namespace SimpleModal;

static Eigen::Vector3d samplePointInTet(
    const Eigen::MatrixXd &TV,
    const Eigen::Vector4i &tet,
    double shrink_factor,
    std::mt19937 &gen)
{
    // Generate random barycentric coordinates
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    double r1 = dist(gen);
    double r2 = dist(gen);
    double r3 = dist(gen);

    // ensure sum of barycentric coordinates is <= 1
    if (r1 + r2 + r3 > 1.0)
    {
        r1 = 1.0 - r1;
        r2 = 1.0 - r2;
        r3 = 1.0 - r3;
    }

    // compute fourth barycentric coordinate
    // such that r1 + r2 + r3 + r4 = 1
    double r4 = 1.0 - r1 - r2 - r3;

    // Sample point using barycentric coordinates
    Eigen::Vector3d p = 
        r1 * TV.row(tet[0]) +
        r2 * TV.row(tet[1]) +
        r3 * TV.row(tet[2]) +
        r4 * TV.row(tet[3]);

    // Shrink p toward origin to avoid sampling too close to the surface
    return shrink_factor * p;
}


// return N x 3 MatrixXd
void Sampling::samplePointsInTetMesh(
    Eigen::MatrixXd &points,
    const TetMesh &mesh,
    int num_points,
    unsigned int seed) 
{
    points.resize(num_points, 3);

    std::mt19937 gen(seed);
    std::uniform_real_distribution<double> dist(0.0, 1.0);

    // Create a scaled copy of mesh.TV
    double global_scale = 0.9;
    Eigen::MatrixXd scaledTV = mesh.TV * global_scale;

    // Compute per-tet volume and cumulative distribution
    std::vector<double> volumes(mesh.TT.rows());
    std::vector<double> cdf(mesh.TT.rows());

    double total_volume = 0.0;
    for (int i = 0; i < mesh.TT.rows(); ++i)
    {
        volumes[i] = FEM::computeTetrahedronVolume(scaledTV, mesh.TT.row(i));
        total_volume += volumes[i];
    }

    // Normalize to CDF
    cdf[0] = volumes[0] / total_volume;
    for (int i = 1; i < mesh.TT.rows(); ++i)
    {
        cdf[i] = cdf[i - 1] + volumes[i] / total_volume;
    }

    // Sample points
    // Eigen::MatrixXd points(num_points, 3);

    for (int i = 0; i < num_points; ++i)
    {
        double r = dist(gen);
        auto it = std::lower_bound(cdf.begin(), cdf.end(), r);
        int tet_index = std::distance(cdf.begin(), it);

        Eigen::Vector3d p = samplePointInTet(scaledTV, mesh.TT.row(tet_index), 0.85, gen);
        points.row(i) = p;
    }
}



void Sampling::samplePointsOnOffsetSurface(
    Eigen::MatrixXd &points,
    const Eigen::MatrixXd &V,
    const Eigen::MatrixXi &F,
    double offset,
    int num_points,
    int offset_surface_resolution,
    unsigned int seed)
{
    int s = offset_surface_resolution;

    Eigen::Vector3d bbox = V.colwise().maxCoeff() - V.colwise().minCoeff();
    double min_side = bbox.minCoeff();
    double offset_val = -offset * min_side;

    // compute offset surface
    Eigen::MatrixXd Vout, GV, S;
    Eigen::MatrixXi Fout;
    Eigen::Vector3i side;

    igl::offset_surface(
        V, F,
        offset_val,
        s,
        igl::SIGNED_DISTANCE_TYPE_PSEUDONORMAL,
        Vout, Fout, GV, side, S);

    // sample points on offset surface
    Eigen::MatrixXd BC;
    Eigen::VectorXi FI;
    std::mt19937 gen(seed);

    igl::random_points_on_mesh(
        num_points,
        Vout, Fout,
        BC, FI, points,
        gen);
}