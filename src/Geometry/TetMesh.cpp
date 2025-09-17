#include "Geometry/TetMesh.h"
#include "Utils/utils.h"
#include <igl/readOBJ.h>
#include <igl/readPLY.h>
#include <igl/readOFF.h>

using namespace SimpleModal;

TetMesh::TetMesh(const Eigen::MatrixXd &TV, const Eigen::MatrixXi &TT, const Eigen::MatrixXi &TF)
    : TV(TV), TT(TT), TF(TF)
{
    buildVertexMapping();
}

TetMesh TetMesh::getTetMeshFromSurfaceMesh(const std::string &surface_mesh_path)
{
    // Read the surface mesh
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    if (surface_mesh_path.find(".obj") != std::string::npos)
    {
        if (!igl::readOBJ(surface_mesh_path, V, F))
        {
            Utils::printError("[TetMesh::getTetrahedralMesh] Failed to read obj file");
            exit(1);
        }
    }
    else if (surface_mesh_path.find(".ply") != std::string::npos)
    {
        if (!igl::readPLY(surface_mesh_path, V, F))
        {
            Utils::printError("[TetMesh::getTetrahedralMesh] Failed to read ply file");
            exit(1);
        }
    }
    else if (surface_mesh_path.find(".off") != std::string::npos)
    {
        if (!igl::readOFF(surface_mesh_path, V, F))
        {
            Utils::printError("[TetMesh::getTetrahedralMesh] Failed to read off file");
            exit(1);
        }
    }
    else
    {
        Utils::printError("[TetMesh::getTetrahedralMesh] Invalid mesh file extension (not obj, ply), off");
        exit(1);
    }

    Eigen::MatrixXd TV; //  vertex positions of the tetrahedral mesh.
    Eigen::MatrixXi TT; //  tetrahedra indices (each row is a tetrahedron defined by 4 vertex indices).
    Eigen::MatrixXi TF; //  face indices.
    // Tetrahedralize the surface mesh
    igl::copyleft::tetgen::tetrahedralize(V, F, "pq1.414Y", TV, TT, TF);

    TetMesh tet_mesh(TV, TT, TF);

    // sanity check V and TV_surface are the same
    for (int i = 0; i < tet_mesh.TV_surface.rows(); i++)
    {
        assert((tet_mesh.TV_surface.row(i) - V.row(i)).norm() < 1e-10);
    }

    return tet_mesh;
}

void TetMesh::buildVertexMapping()
{
    // Get the surface vertices
    std::set<int> surface_vertices_idx;
    for (int i = 0; i < TF.rows(); i++)
    {
        surface_vertices_idx.insert(TF(i, 0));
        surface_vertices_idx.insert(TF(i, 1));
        surface_vertices_idx.insert(TF(i, 2));
    }
    TV_surface = Eigen::MatrixXd::Zero(surface_vertices_idx.size(), 3);

    int surface_vert_count = 0;
    for (int i = 0; i < TV.rows(); i++)
    {
        if (surface_vertices_idx.count(i))
        {
            _tet_vert_to_surface_vert_map[i] = surface_vert_count;
            _surface_vert_to_tet_vert_map[surface_vert_count] = i;
            TV_surface.row(surface_vert_count) = TV.row(i);
            surface_vert_count++;
        }
    }
    std::cout << "_tet_vert_to_surface_vert_map.size() = " << _tet_vert_to_surface_vert_map.size() << std::endl;
    std::cout << "_surface_vert_to_tet_vert_map.size() = " << _surface_vert_to_tet_vert_map.size() << std::endl;
}

int TetMesh::mapSurface2Tet(int surface_vertex_idx) const
{
    if (_surface_vert_to_tet_vert_map.find(surface_vertex_idx) == _surface_vert_to_tet_vert_map.end())
    {
        Utils::printError("Surface vertex index not found in the map");
        exit(1);
    }
    return _surface_vert_to_tet_vert_map.at(surface_vertex_idx);
}
int TetMesh::mapTet2Surface(int tet_vertex_idx) const
{
    if (_tet_vert_to_surface_vert_map.find(tet_vertex_idx) == _tet_vert_to_surface_vert_map.end())
    {
        Utils::printError("Tet vertex index not found in the map");
        exit(1);
    }
    return _tet_vert_to_surface_vert_map.at(tet_vertex_idx);
}