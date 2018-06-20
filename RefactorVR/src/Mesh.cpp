#include "Mesh.h"

Mesh::Mesh(Eigen::MatrixXd & V_, Eigen::MatrixXi & F_, Eigen::VectorXi & vertex_boundary_markers_, Eigen::VectorXi &part_of_original_stroke_, Eigen::VectorXi &new_mapped_indices_, Eigen::VectorXi &sharp_edge_, int ID_): V(V_), F(F_), vertex_boundary_markers(vertex_boundary_markers_), part_of_original_stroke(part_of_original_stroke_), new_mapped_indices(new_mapped_indices_), sharp_edge(sharp_edge_), ID(ID_)  {
}

Mesh::Mesh() :
	V(Eigen::MatrixXd(0, 3)),
	F(Eigen::MatrixXi(0, 3)),
	vertex_boundary_markers(Eigen::VectorXi()),
	part_of_original_stroke(Eigen::VectorXi()),
	new_mapped_indices(Eigen::VectorXi()),
	sharp_edge(Eigen::VectorXi()),
	ID(-1),
	patches(),
	face_patch_map()
	{
	/*V = (Eigen::MatrixXd(0, 3));
	F = (Eigen::MatrixXi(0, 3));
	vertex_boundary_markers = (Eigen::VectorXi());
	part_of_original_stroke = (Eigen::VectorXi());
	new_mapped_indices = (Eigen::VectorXi());
	sharp_edge = (Eigen::VectorXi());
	ID = (-1);*/
	}

Mesh & Mesh::operator=(Mesh other){
	/*std::swap(V, other.V);
	std::swap(F, other.F);
	std::swap(vertex_boundary_markers, other.vertex_boundary_markers);
	std::swap(part_of_original_stroke, other.part_of_original_stroke);
	std::swap(new_mapped_indices, other.new_mapped_indices);
	std::swap(sharp_edge, other.sharp_edge);
	std::swap(ID, other.ID);*/
	F.resize(other.F.rows(), 3);
	F = other.F;
	vertex_boundary_markers = other.vertex_boundary_markers;
	part_of_original_stroke = other.part_of_original_stroke;
	new_mapped_indices = other.new_mapped_indices;
	sharp_edge = other.sharp_edge;
	ID = other.ID;
	patches = other.patches;
	face_patch_map = other.face_patch_map;
	return *this;
}

Mesh::~Mesh() {
}