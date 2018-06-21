#include "Patch.h"
#include <igl/edge_topology.h>

static Eigen::MatrixXi EV, FE, EF;
static int next_mesh_ID;

/*Patch& Patch::operator=(Patch other) {
	mesh = other.mesh;
	parent_vertices = other.parent_vertices;
	parent_faces = other.parent_faces;
	new_vertex_boundary_markers = other.new_vertex_boundary_markers;
	new_part_of_original_stroke = other.new_part_of_original_stroke;
	//new_mapped_indices = other.new_mapped_indices;
	return *this;
}*/

std::vector<Patch*> Patch::init_patches(Mesh& h) {
	next_mesh_ID = h.ID + 1;
	igl::edge_topology(h.V, h.F, EV, FE, EF);

	std::vector<Patch*> patches;
	h.face_patch_map.resize(h.F.rows(), nullptr);
	h.mesh_to_patch_indices.resize(h.V.rows());
	h.mesh_to_patch_indices.setConstant(-1);
	
	for (int i = 0; i < h.F.rows(); i++) {
		Patch* face_patch = h.face_patch_map[i];
		if (face_patch == nullptr) {
			Patch* new_patch = new Patch();
			Eigen::VectorXi faces;
			propagate_patch(new_patch, i, faces, h.face_patch_map, h.sharp_edge);
			(*new_patch).create_mesh_structure(h, faces);
			patches.push_back(new_patch);
		}
	}
	return patches;
}

void Patch::propagate_patch(Patch* patch, int face, Eigen::VectorXi& faces, std::vector<Patch*> &face_patch_map, Eigen::VectorXi& sharp_edge) {
	face_patch_map[face] = patch;
	faces.conservativeResize(faces.rows() + 1);
	faces.tail(1) << face;
	for (int i = 0; i < 3; i++) {
		int edge = FE(face, i);
		int face2 = (EF(edge, 0) == face) ? EF(edge, 1) : EF(edge, 0); //get the adjacent polygon
		if (face_patch_map[face2] == nullptr && !sharp_edge[edge]) {
			propagate_patch(patch, face2, faces, face_patch_map, sharp_edge);
		}
	}
}

void Patch::create_mesh_structure(Mesh& m, Eigen::VectorXi& faces) {
	patch_vertices.resize(0, 3);
	patch_vertex_is_init = Eigen::VectorXi::Zero(m.V.rows());
	patch_faces.resize(0,3);
	patch_edge_is_init = Eigen::VectorXi::Zero(EF.rows());
	mesh_to_patch_indices = Eigen::VectorXi::Constant(m.V.rows(),-1);

	for (int i = 0; i < faces.rows(); i++) {
		for (int j = 0; j < 3; j++) {
			int edge = FE(faces[i], j);
			get_patch_edge(edge, patch_edge_is_init, patch_vertex_is_init, faces[i], patch_vertices, m.sharp_edge, new_sharp_edge, m.V, m.vertex_boundary_markers, m.part_of_original_stroke, m.new_mapped_indices, mesh_to_patch_indices);
		}
	
		parent_faces.conservativeResize(parent_faces.rows() + 1, Eigen::NoChange);
		parent_faces[i] = faces[i];	
	}

	for (int i = 0; i < faces.rows(); i++) {
		patch_faces.conservativeResize(patch_faces.rows() + 1, Eigen::NoChange);
		patch_faces.bottomRows(1) << mesh_to_patch_indices[m.F(faces[i], 0)], mesh_to_patch_indices[m.F(faces[i], 1)], mesh_to_patch_indices[m.F(faces[i], 2)];
	}
	mesh.V = patch_vertices;
	mesh.F = patch_faces;
	mesh.vertex_boundary_markers = new_vertex_boundary_markers;
	mesh.part_of_original_stroke = new_part_of_original_stroke;
	mesh.new_mapped_indices = new_mapped_indices;
	mesh.sharp_edge = new_sharp_edge;
	mesh.ID = next_mesh_ID;
	mesh.mesh_to_patch_indices = mesh_to_patch_indices; //Maps a vertex index from the base mesh to a vertex index in the mesh of a underlying patch
	next_mesh_ID++;
}

void Patch::get_patch_edge(int edge, Eigen::VectorXi& patch_edge_is_init, Eigen::VectorXi& patch_vertex_is_init, int face, Eigen::MatrixXd& patch_vertices, Eigen::VectorXi& sharp_edge, Eigen::VectorXi& new_sharp_edge, Eigen::MatrixXd& V_orig, Eigen::VectorXi& boundary_markers_orig, Eigen::VectorXi& part_of_original_orig, Eigen::VectorXi& new_mapped_indices_orig, Eigen::VectorXi& mesh_to_patch_indices) {
	if (sharp_edge(edge)) {
		//Do nothing
	}
	else if (patch_edge_is_init[edge]) {
		return;
	}

	int start = EV(edge, 0);
	int end = EV(edge, 1);
	get_patch_vertex(start, face, patch_vertices, patch_vertex_is_init, V_orig, boundary_markers_orig, part_of_original_orig, new_mapped_indices_orig, mesh_to_patch_indices);
	get_patch_vertex(end, face, patch_vertices, patch_vertex_is_init, V_orig, boundary_markers_orig, part_of_original_orig, new_mapped_indices_orig, mesh_to_patch_indices);
	new_sharp_edge.conservativeResize(new_sharp_edge.rows() + 1);
	new_sharp_edge.tail(1) << sharp_edge[edge];
	patch_edge_is_init[edge] = 1;

}


void Patch::get_patch_vertex(int v_idx, int face, Eigen::MatrixXd& patch_vertices, Eigen::VectorXi& patch_vertex_is_init, Eigen::MatrixXd& V_orig, Eigen::VectorXi& boundary_markers_orig, Eigen::VectorXi& part_of_original_orig, Eigen::VectorXi& new_mapped_indices_orig, Eigen::VectorXi& mesh_to_patch_indices) {
	if (patch_vertex_is_init[v_idx]) {
		return;
	}

	patch_vertices.conservativeResize(patch_vertices.rows() + 1, Eigen::NoChange);
	patch_vertices.bottomRows(1) = V_orig.row(v_idx);
	parent_vertices.conservativeResize(parent_vertices.rows() + 1);
	parent_vertices.tail(1) << v_idx;
	mesh_to_patch_indices[v_idx] = parent_vertices.rows() - 1;
	new_vertex_boundary_markers.conservativeResize(new_vertex_boundary_markers.rows() + 1);
	new_vertex_boundary_markers.tail(1) << boundary_markers_orig.row(v_idx);
	new_part_of_original_stroke.conservativeResize(new_part_of_original_stroke.rows() + 1);
	new_part_of_original_stroke.tail(1) << part_of_original_orig.row(v_idx);
	new_mapped_indices.conservativeResize(new_mapped_indices.rows() + 1);
	new_mapped_indices.tail(1) << -1;// new_mapped_indices_orig.row(v_idx); //TODO: CHECK IF WE NEED THE NEW_MAPPED_INDICES OR IF THESE ARE ALWAYS UNINITIALIZED WHEN THIS IS CALLED
	patch_vertex_is_init[v_idx] = 1;
}

void Patch::update_parent_vertex_positions(Eigen::MatrixXd& base_V) {
	std::cout << " Base size: " << base_V.rows() << std::endl;
	std::cout << "other size: " << mesh.V.rows() << std::endl;
	for (int i = 0; i < mesh.V.rows(); i++) {
		base_V.row(parent_vertices[i]) = mesh.V.row(i);
	}
}