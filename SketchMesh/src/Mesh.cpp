#include "Mesh.h"

Mesh::Mesh(Eigen::MatrixXd & V_, Eigen::MatrixXi & F_, const Eigen::VectorXi & vertex_boundary_markers_, const Eigen::VectorXi &part_of_original_stroke_, int ID_): V(V_), F(F_), vertex_boundary_markers(vertex_boundary_markers_), part_of_original_stroke(part_of_original_stroke_), ID(ID_)  {

}

Mesh::~Mesh() {
}

/*bool Mesh::operator==(const Mesh &other) const {
	return (V == other.V
		&& F == other.F
		&& vertex_boundary_markers == other.vertex_boundary_markers);
}*/
