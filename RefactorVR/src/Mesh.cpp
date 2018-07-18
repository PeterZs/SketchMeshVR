#include "Mesh.h"

Mesh::Mesh(Eigen::MatrixXd & V_, Eigen::MatrixXi & F_, Eigen::VectorXi & vertex_boundary_markers_, Eigen::VectorXi & edge_boundary_markers_, Eigen::VectorXi& vertex_is_fixed_, Eigen::VectorXi &part_of_original_stroke_, Eigen::VectorXi &new_mapped_indices_, Eigen::VectorXi &sharp_edge_, int ID_): V(V_), F(F_), vertex_boundary_markers(vertex_boundary_markers_), edge_boundary_markers(edge_boundary_markers_), vertex_is_fixed(vertex_is_fixed_), part_of_original_stroke(part_of_original_stroke_), new_mapped_indices(new_mapped_indices_), sharp_edge(sharp_edge_), ID(ID_)  {
}

Mesh::Mesh() :
	V(Eigen::MatrixXd(0, 3)),
	F(Eigen::MatrixXi(0, 3)),
	vertex_boundary_markers(Eigen::VectorXi()),
	edge_boundary_markers(Eigen::VectorXi()),
	vertex_is_fixed(Eigen::VectorXi()),
	part_of_original_stroke(Eigen::VectorXi()),
	new_mapped_indices(Eigen::VectorXi()),
	sharp_edge(Eigen::VectorXi()),
	ID(-1),
	patches(),
	face_patch_map(),
	mesh_to_patch_indices()
	{}

Mesh & Mesh::operator=(Mesh other){
	F.resize(other.F.rows(), 3);
	F = other.F;
	vertex_boundary_markers = other.vertex_boundary_markers;
	edge_boundary_markers = other.edge_boundary_markers;
	vertex_is_fixed = other.vertex_is_fixed;
	part_of_original_stroke = other.part_of_original_stroke;
	new_mapped_indices = other.new_mapped_indices;
	sharp_edge = other.sharp_edge;
	ID = other.ID;
	patches = other.patches;
	face_patch_map = other.face_patch_map;
	mesh_to_patch_indices = other.mesh_to_patch_indices;
	return *this;
}

Mesh::~Mesh() {
}