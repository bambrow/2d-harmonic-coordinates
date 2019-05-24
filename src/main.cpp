// Geometric Modeling
// Final Project
// 2D Harmonic Coordinates
// Author: Weiqiang Li
// wl1731@nyu.edu

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/triangle/triangulate.h>
#include <igl/project.h>
#include <igl/unproject.h>
#include <igl/readOFF.h>
#include <igl/slice.h>
#include <igl/cotmatrix.h>
#include <igl/boundary_facets.h>
#include <igl/unique.h>
#include <string>
#include <fstream>
#include <sstream>
#include <unordered_map>

// this project is aimed for harmonic coordinates for mesh deformation
// the supposed usage is any of the following:
// ./final_bin
// which will run the algorithm based on default mesh and cage
// ./final_bin <mesh.off>
// which will read the mesh and automatically build a cage based on the boundary of the mesh
// ./final_bin <mesh.off> <cage.cage>
// which will read the mesh and cage from file
// where mesh.off is the mesh information (V and F)
// in mesh.off, the z coordinates are assumed 0 (not used)
// cage.cage is the cage information (V)
// in cage.cage, the z coordinates are assumed 0 (not used)
// the vertices in cage.cage should be either ordered clockwise or counterclockwise
// interaction manual:
// use mouse (click and drag) to change the location of cage vertex
// the new mesh will be calculated automatically
// press R to reset the cage and mesh to original position
// press U to undo the cage vertex change (only once)
// use W/A/S/D to control the cage vertex location after selected by mouse
// the new mesh will be cauculated automatically

using namespace std;
using namespace Eigen;
using Viewer = igl::opengl::glfw::Viewer;

// debug flags
// #define DEBUG_1
// #define DEBUG_2
// #define DEBUG_3

// mesh
MatrixXd V;
MatrixXi F;
// mesh archive
MatrixXd V1;
MatrixXd V2;

// cage
MatrixXd CV;
MatrixXi CE;
// cage archive
MatrixXd CV1;
MatrixXd CV2;

// triangulate helper
MatrixXd VV0;
MatrixXi VF0;

// harmonic helper
SparseMatrix<double> Aff, Afc;
MatrixXd H;

// interactions helper
bool mouse_down = false;
int current_cage_vertex = -1;
double keyboard_stride = 0.0;

// declarations
void plot_mesh_and_cage(igl::opengl::glfw::Viewer &viewer);
bool callback_key_down(Viewer &viewer, unsigned char key, int modifiers);
int find_nearest_cage(RowVector3d loc);
bool callback_mouse_down(Viewer& viewer, int, int);
bool callback_mouse_move(Viewer& viewer, int mouse_x, int mouse_y);
bool callback_mouse_up(Viewer& viewer, int, int);
void solve_prepare();
void solve_harmonic();
bool read_cage_from_file(string cage_filename);
void generate_cage();


// plots mesh (in V and F) and cage (in CV and CE)
// during mouse click and drag, the selected cage vertex will be marked green
void plot_mesh_and_cage(igl::opengl::glfw::Viewer &viewer) {
  viewer.data().clear();
  viewer.data().set_mesh(V, F);

  if (mouse_down) {
    for (unsigned i = 0; i < CV.rows(); i++) {
      if (i == current_cage_vertex) {
        viewer.data().add_points(CV.row(i), RowVector3d(0,1,0));
      } else {
        viewer.data().add_points(CV.row(i), RowVector3d(1,0,0));
      }
    }
  } else {
    viewer.data().add_points(CV, RowVector3d(1,0,0));
  }
  for (unsigned i = 0; i < CE.rows(); i++) {
    viewer.data().add_edges(
      CV.row(CE(i,0)),
      CV.row(CE(i,1)),
      RowVector3d(1,0,0)
    );
  }
}

// key_down callback
// accepts R (reset), U (undo), W/A/S/D (cage vertex move)
bool callback_key_down(Viewer &viewer, unsigned char key, int modifiers) {
  if (key == 'R') {
    CV = CV1;
    V = V1;
    plot_mesh_and_cage(viewer);
    return true;
  }
  else if (key == 'U') {
    CV = CV2;
    V = V2;
    plot_mesh_and_cage(viewer);
    return true;
  }
  else if (key == 'W' || key == 'S' || key == 'A' || key == 'D') {
    if (current_cage_vertex != -1) {
      if (key == 'W') {
        CV(current_cage_vertex, 1) += keyboard_stride;
      } 
      else if (key == 'S') {
        CV(current_cage_vertex, 1) -= keyboard_stride;
      }
      else if (key == 'A') {
        CV(current_cage_vertex, 0) -= keyboard_stride;
      }
      else {
        CV(current_cage_vertex, 0) += keyboard_stride;
      }
      CV2 = CV;
      V2 = V;
      solve_harmonic();
      #ifdef DEBUG_3
      cout << V << endl << endl;
      #endif
      plot_mesh_and_cage(viewer);
    }
    return true;
  }
  return false;
}

// find the nearest cage vertex based on the unprojected mouse location
// the location should be close to at least one of the cage vertices
// otherwise will return -1
int find_nearest_cage(RowVector3d loc) {
  int nearest = -1;
  double min_distance = numeric_limits<double>::max();
  for (unsigned i = 0; i < CV.rows(); i++) {
    RowVector3d diff = loc - CV.row(i);
    if (diff.norm() < min_distance) {
      min_distance = diff.norm();
      nearest = i;
    }
  }
  RowVector3d cage_max = CV.colwise().maxCoeff();
  RowVector3d cage_min = CV.colwise().minCoeff();
  #ifdef DEBUG_1
  cout << cage_max << endl;
  cout << cage_min << endl;
  #endif
  double distance_threshold = (cage_max(0) - cage_min(0) + cage_max(1) - cage_min(1)) / 2 * 0.1;
  keyboard_stride = distance_threshold / 2;
  if (min_distance > distance_threshold) {
    nearest = -1;
  }
  return nearest;
}

// mouse_down callback
// unproject the mouse location and find the a cage vertex nearby
// and then prepare for the drag
bool callback_mouse_down(Viewer& viewer, int, int) {
  int mx = viewer.current_mouse_x;
  int my = viewer.core.viewport(3) - viewer.current_mouse_y;
  RowVector3d origin_point;
  origin_point.setZero();
  RowVector3d projected;
  igl::project(origin_point, viewer.core.view, viewer.core.proj, viewer.core.viewport, projected);
  #ifdef DEBUG_3
  cout << "projected: " << projected << endl;
  #endif
  double mz = projected[2];
  RowVector3d target;
  target.setZero();
  igl::unproject(RowVector3d(mx,my,mz), viewer.core.view, viewer.core.proj, viewer.core.viewport, target);
  int cid = find_nearest_cage(target);
  #ifdef DEBUG_3
  cout << "cage selected: " << cid << endl;
  cout << "mouse location: " << mx << ", " << my << "," << mz << endl;
  cout << "target location: " << target << endl;
  #endif
  if (cid == -1) {
    return false;
  }
  current_cage_vertex = cid;
  mouse_down = true;
  CV2 = CV;
  V2 = V;
  plot_mesh_and_cage(viewer);
  return true;
}

// mouse_move callback
// works only when mouse is down (when dragging)
// unproject the mouse location and move the selected cage vertex accordingly
// and then solve the system 
bool callback_mouse_move(Viewer& viewer, int mouse_x, int mouse_y) {
  if (mouse_down) {
    int mx = mouse_x;
    int my = viewer.core.viewport(3) - mouse_y;
    RowVector3d origin_point;
    origin_point.setZero();
    RowVector3d projected;
    igl::project(origin_point, viewer.core.view, viewer.core.proj, viewer.core.viewport, projected);
    #ifdef DEBUG_3
    cout << "projected: " << projected << endl;
    #endif
    double mz = projected[2];
    RowVector3d target;
    target.setZero();
    igl::unproject(RowVector3d(mx,my,mz), viewer.core.view, viewer.core.proj, viewer.core.viewport, target);
    #ifdef DEBUG_3
    cout << "mouse location: " << mx << ", " << my << "," << mz << endl;
    cout << "target location: " << target << endl;
    #endif
    CV.row(current_cage_vertex) = target;
    solve_harmonic();
    plot_mesh_and_cage(viewer);
  }
  return true;
}

// mouse_up callback
bool callback_mouse_up(Viewer& viewer, int, int) {
  mouse_down = false;
  plot_mesh_and_cage(viewer);
  return true;
}

// solver prepare
// use Laplace equation to solve the system
// stores the harmonic matrix in H
// consider mesh and cage as a whole new mesh
// and use cage vertices as boundary vertices
// then utilize variable elimination to build the linear system
void solve_prepare() {
  MatrixXd VV(CV.rows()+V.rows(), 2);
  for (unsigned i = 0; i < CV.rows(); i++) {
    VV(i,0) = CV(i,0);
    VV(i,1) = CV(i,1);
  }
  for (unsigned i = 0; i < V.rows(); i++) {
    VV(CV.rows()+i,0) = V(i,0);
    VV(CV.rows()+i,1) = V(i,1);
  }
  MatrixXd H0(0, 2);
  igl::triangle::triangulate(VV,CE,H0,"Q",VV0,VF0);
  VV0.conservativeResize(VV0.rows(),3);
  VV0.col(2).setZero();
  #ifdef DEBUG_1
  cout << "VV = " << endl << VV << endl;
  cout << "VV0 = " << endl << VV0 << endl;
  #endif
  SparseMatrix<double> L;
  igl::cotmatrix(VV0,VF0,L);
  VectorXi all, in, b;
  igl::colon<int>(0,VV0.rows()-1,all);
  igl::colon<int>(CV.rows(),VV0.rows()-1,in);
  igl::colon<int>(0,CV.rows()-1,b);
  SparseMatrix<double> A = L * (-1);
  // SparseMatrix<double> Aff, Afc;
  igl::slice(A,in,in,Aff);
  igl::slice(A,in,b,Afc);
  
  SimplicialLLT<SparseMatrix<double>> solver;
  solver.compute(Aff);
  H = solver.solve(MatrixXd(Afc)*(-1));
  
  #ifdef DEBUG_2
  cout << "VV0: " << VV0.rows() << " * " << VV0.cols() << endl;
  cout << "VF0: " << VF0.rows() << " * " << VF0.cols() << endl;
  cout << "L: " << L.rows() << " * " << L.cols() << endl;
  cout << "A: " << A.rows() << " * " << A.cols() << endl;
  cout << "all: " << all.rows() << " * " << all.cols() << endl;
  cout << "in: " << in.rows() << " * " << in.cols() << endl;
  cout << "b: " << b.rows() << " * " << b.cols() << endl;
  cout << "Aff: " << Aff.rows() << " * " << Aff.cols() << endl;
  cout << "Afc: " << Afc.rows() << " * " << Afc.cols() << endl;
  cout << "H: " << H.rows() << " * " << H.cols() << endl;
  #endif
}

// compute the new mesh based on H
void solve_harmonic() {
  /*
  for (unsigned i = 0; i < V.cols(); i++) {
    VectorXd bc = CV.col(i);
    SimplicialLLT<SparseMatrix<double>> solver(Aff);
    VectorXd XX = solver.solve(MatrixXd(Afc) * (-1) * bc);
    V.col(i) = XX;
  }
  */
  
  MatrixXd NV = H * CV;
  for (unsigned i = 0; i < NV.rows(); i++) {
    V.row(i) = NV.row(i);
  }
  
}

// read cage information from .cage file
// if failed, the main function will use automatic cage generation instead
bool read_cage_from_file(string cage_filename) {
  ifstream in(cage_filename);
  if (!in.is_open()) {
    cout << "Error: cannot open the cage file! Automatic cage generation will be used." << endl;
    return false;
  }
  CV.resize(8, 3);
  int index = 0;
  string line;
  while (getline(in, line)) {
    index++;
    istringstream iss(line);
    if (index > CV.rows()) {
      CV.conservativeResize(CV.rows()*2, CV.cols());
    }
    for (unsigned i = 0; i < 3; i++) {
      double xyz;
      if (iss >> xyz) {
        CV(index-1,i) = xyz; 
      } else {
        cout << "Error: cannot open the cage file! Automatic cage generation will be used." << endl;
        return false;
      }
    }
  }
  CV.conservativeResize(index, CV.cols());
  in.close();
  return true;
}

// automatically generate cage based on mesh boundary
// compute the mesh centroid and move the boundary vertices even further to build the cage
void generate_cage() {
  int vn = V.rows();
  RowVector3d mesh_centriod;
  mesh_centriod.setZero();
  for (unsigned i = 0; i < vn; i++) {
    mesh_centriod += V.row(i);
  }
  mesh_centriod /= vn;
  MatrixXi VE;
  igl::boundary_facets(F,VE);
  #ifdef DEBUG_2
  cout << "VE: " << endl << VE << endl << endl;
  #endif
  VectorXi CC, IA, IC;
  igl::unique(VE,CC,IA,IC);
  #ifdef DEBUG_2
  cout << "CC: " << endl << CC << endl << endl;
  cout << "IA: " << endl << IA << endl << endl;
  cout << "IC: " << endl << IC << endl << endl;
  #endif
  MatrixXd VI;
  igl::slice(V,CC,1,VI);
  CV.resizeLike(VI);
  for (unsigned i = 0; i < VI.rows(); i++) {
    CV.row(i) = (VI.row(i) - mesh_centriod) * 1.5 + mesh_centriod;
  }
  unordered_map<int,int> dict;
  for (unsigned i = 0; i < CC.rows(); i++) {
    dict[CC(i)] = i;
  }
  CE.resizeLike(VE);
  for (unsigned i = 0; i < VE.rows(); i++) {
    CE(i,0) = dict[VE(i,0)]; 
    CE(i,1) = dict[VE(i,1)];
  }
}

// main function
int main(int argc, char *argv[])
{
  if (argc <= 1) {

    V.resize(4,3);
    V << 0,0,0,
         1,0,0,
         0,1,0,
         1,1,0;
    F.resize(2,3);
    F << 0,1,2,
         1,3,2;
    CV.resize(4,3);
    CV << -0.5,-0.5,0,
          1.5,-0.5,0,
          1.5,1.5,0,
          -0.5,1.5,0;
    CE.resize(4,2);
    CE << 0,1,
          1,2,
          2,3,
          3,0;

  }
  else if (argc == 2) {

    string mesh_file_name = argv[1];
    igl::readOFF(mesh_file_name,V,F);
    generate_cage();

  }
  else if (argc >= 3) {

    string mesh_file_name = argv[1];
    string mesh_cage_name = argv[2];

    igl::readOFF(mesh_file_name,V,F);
    if (!read_cage_from_file(mesh_cage_name)) {
      generate_cage();
    } else {
      CE.resize(CV.rows(),2);
      for (unsigned i = 0; i < CV.rows(); i++) {
        CE(i,0) = i;
        CE(i,1) = (i+1)%CV.rows();
      }
    }

  }
  else {
    cout << "Usage: ./final_bin <mesh.off> <cage.cage>" << endl;
    cout << "Automatic cage generation: " << endl;
    cout << "Usage: ./final_bin <mesh.off>" << endl;
    cout << "Default mesh and cage: " << endl;
    cout << "Usage: ./final_bin" << endl;
    exit(1);
  }

  CV1 = CV2 = CV;
  V1 = V2 = V;

  // Plot the mesh
  igl::opengl::glfw::Viewer viewer;

  viewer.callback_key_down = &callback_key_down;
  viewer.callback_mouse_down = &callback_mouse_down;
  viewer.callback_mouse_move = &callback_mouse_move;
  viewer.callback_mouse_up = &callback_mouse_up;

  solve_prepare();

  plot_mesh_and_cage(viewer);
  viewer.core.align_camera_center(V,F);

  viewer.launch();

}
