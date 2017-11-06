#include "Mesh.h"

Mesh::Mesh(Eigen::MatrixXd & V_, Eigen::MatrixXi & F_, const Eigen::VectorXi & vertex_boundary_markers_): V(V_), F(F_), vertex_boundary_markers(vertex_boundary_markers_)  {

}

Mesh::~Mesh() {
}
