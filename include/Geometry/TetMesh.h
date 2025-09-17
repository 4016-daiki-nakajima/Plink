#pragma once

#include <Eigen/Dense>
#include <string>
#include <igl/copyleft/tetgen/tetrahedralize.h>
#include <map>

namespace SimpleModal
{
    class TetMesh
    {
    public:
        TetMesh() {};
        TetMesh(const Eigen::MatrixXd& TV, const Eigen::MatrixXi& TT, const Eigen::MatrixXi& TF);

        static TetMesh getTetMeshFromSurfaceMesh(const std::string& surface_mesh_path);

        // Tetrahedralized interior
        Eigen::MatrixXd TV;
        Eigen::MatrixXi TT;
        Eigen::MatrixXi TF;
        Eigen::MatrixXd TV_surface;

        int mapSurface2Tet(int surface_vertex_idx) const;
        int mapTet2Surface(int tet_vertex_idx) const;

        void buildVertexMapping(); 

    protected:
        std::map<int, int> _tet_vert_to_surface_vert_map; 
        std::map<int, int> _surface_vert_to_tet_vert_map; 
    };
};
