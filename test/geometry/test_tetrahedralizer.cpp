#include "TetMesh.h"
#include <igl/barycenter.h>
#include <igl/opengl/glfw/Viewer.h>
#include <iostream>

// Ref: https://github.com/libigl/libigl/blob/main/tutorial/605_Tetgen/main.cpp
Eigen::MatrixXd TV;
Eigen::MatrixXi TT;
Eigen::MatrixXi TF;
Eigen::MatrixXd B;

using namespace SimpleModal;

// This function is called every time a keyboard button is pressed
bool key_down(igl::opengl::glfw::Viewer &viewer, unsigned char key,
              int modifier) {
  using namespace std;
  using namespace Eigen;

  if (key >= '1' && key <= '9') {
    double t = double((key - '1') + 1) / 9.0;

    VectorXd v = B.col(2).array() - B.col(2).minCoeff();
    v /= v.col(0).maxCoeff();

    vector<int> s;

    for (unsigned i = 0; i < v.size(); ++i)
      if (v(i) < t)
        s.push_back(i);

    MatrixXd V_temp(s.size() * 4, 3);
    MatrixXi F_temp(s.size() * 4, 3);

    for (unsigned i = 0; i < s.size(); ++i) {
      V_temp.row(i * 4 + 0) = TV.row(TT(s[i], 0));
      V_temp.row(i * 4 + 1) = TV.row(TT(s[i], 1));
      V_temp.row(i * 4 + 2) = TV.row(TT(s[i], 2));
      V_temp.row(i * 4 + 3) = TV.row(TT(s[i], 3));
      F_temp.row(i * 4 + 0) << (i * 4) + 0, (i * 4) + 1, (i * 4) + 3;
      F_temp.row(i * 4 + 1) << (i * 4) + 0, (i * 4) + 2, (i * 4) + 1;
      F_temp.row(i * 4 + 2) << (i * 4) + 3, (i * 4) + 2, (i * 4) + 0;
      F_temp.row(i * 4 + 3) << (i * 4) + 1, (i * 4) + 2, (i * 4) + 3;
    }

    viewer.data().clear();
    viewer.data().set_mesh(V_temp, F_temp);
    viewer.data().set_face_based(true);
  }

  return false;
}
int main() {
  const std::string root_path = "C:/Users/Zhehao/Research/SimpleModal/";
  TetMesh tet_mesh =
      TetMesh::getTetMeshFromSurfaceMesh(root_path + "asset/fertility.off");
  igl::barycenter(tet_mesh.TV, tet_mesh.TT, B);

  TV = tet_mesh.TV;
  TT = tet_mesh.TT;
  TF = tet_mesh.TF;

  igl::opengl::glfw::Viewer viewer;
  viewer.callback_key_down = &key_down;
  key_down(viewer, '5', 0);
  viewer.launch();
  return 0;
}
