#include "Mesh.h"

Mesh::Mesh(Eigen::MatrixXd & V_, Eigen::MatrixXi & F_, Eigen::VectorXi & vertex_boundary_markers_, Eigen::VectorXi &part_of_original_stroke_, int ID_): V(V_), F(F_), vertex_boundary_markers(vertex_boundary_markers_), part_of_original_stroke(part_of_original_stroke_), ID(ID_)  {

}

Mesh::~Mesh() {
}