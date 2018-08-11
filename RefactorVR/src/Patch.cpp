#include "Patch.h"
#include <igl/edge_topology.h>

static Eigen::MatrixXi EV, FE, EF, patchEV, patchFE, patchEF;
static int next_mesh_ID;
static Eigen::MatrixXd Vert;

std::vector<Patch*> Patch::init_patches(Mesh& h) {
	next_mesh_ID = h.ID + 1;
	igl::edge_topology(h.V, h.F, EV, FE, EF);
	Vert = h.V;

	std::vector<Patch*> patches;
	h.face_patch_map.clear();
	h.face_patch_map.resize(h.F.rows(), nullptr);

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
	patch_faces.resize(0, 3);
	patch_edge_is_init = Eigen::VectorXi::Zero(EF.rows());
	mesh_to_patch_indices = Eigen::VectorXi::Constant(m.V.rows(), -1);
	Eigen::MatrixXi tmp_sharp_edge(0, 2), tmp_edge_boundary_markers(0, 3);

	for (int i = 0; i < faces.rows(); i++) {
		for (int j = 0; j < 3; j++) {
			int edge = FE(faces[i], j);
			get_patch_edge(edge, patch_edge_is_init, patch_vertex_is_init, faces[i], patch_vertices, m.sharp_edge, tmp_sharp_edge, m.edge_boundary_markers, tmp_edge_boundary_markers, m.V, m.vertex_is_fixed, m.new_mapped_indices, mesh_to_patch_indices);
		}

		//parent_faces.conservativeResize(parent_faces.rows() + 1, Eigen::NoChange);
	//	parent_faces[i] = faces[i];
	}

	for (int i = 0; i < faces.rows(); i++) {
		patch_faces.conservativeResize(patch_faces.rows() + 1, Eigen::NoChange);
		patch_faces.bottomRows(1) << mesh_to_patch_indices[m.F(faces[i], 0)], mesh_to_patch_indices[m.F(faces[i], 1)], mesh_to_patch_indices[m.F(faces[i], 2)];
	}

	igl::edge_topology(patch_vertices, patch_faces, patchEV, patchFE, patchEF);

	Eigen::VectorXi col1Equals, col2Equals;
	int start, end, val, equal_pos;
	new_sharp_edge.resize(patchEV.rows());
	new_sharp_edge.setZero();
	for (int i = 0; i < tmp_sharp_edge.rows(); i++) { //Go over all sharp edges in the patch and set them
		start = mesh_to_patch_indices[tmp_sharp_edge(i, 0)];
		end = mesh_to_patch_indices[tmp_sharp_edge(i, 1)];
		col1Equals = patchEV.col(0).cwiseEqual(std::min(start, end)).cast<int>();
		col2Equals = patchEV.col(1).cwiseEqual(std::max(start, end)).cast<int>();
		val = (col1Equals + col2Equals).maxCoeff(&equal_pos);
		if (val == 2) {
			new_sharp_edge[equal_pos] = 1;
		}
	}

	new_edge_boundary_markers.resize(patchEV.rows());
	for (int i = 0; i < tmp_edge_boundary_markers.rows(); i++) {
		start = mesh_to_patch_indices[tmp_edge_boundary_markers(i, 0)];
		end = mesh_to_patch_indices[tmp_edge_boundary_markers(i, 1)];
		col1Equals = patchEV.col(0).cwiseEqual(std::min(start, end)).cast<int>();
		col2Equals = patchEV.col(1).cwiseEqual(std::max(start, end)).cast<int>();
		val = (col1Equals + col2Equals).maxCoeff(&equal_pos);
		if (val == 2) {
			new_edge_boundary_markers[equal_pos] = tmp_edge_boundary_markers(i, 2);
		}
	}

	mesh.V = patch_vertices;
	mesh.F = patch_faces;
	mesh.edge_boundary_markers = new_edge_boundary_markers;
	mesh.vertex_is_fixed = new_vertex_is_fixed;
	mesh.new_mapped_indices = new_mapped_indices;
	mesh.sharp_edge = new_sharp_edge;
	mesh.ID = next_mesh_ID;
	mesh.mesh_to_patch_indices = mesh_to_patch_indices; //Maps a vertex index from the base mesh to a vertex index in the mesh of a underlying patch
	next_mesh_ID++;
}

void Patch::get_patch_edge(int edge, Eigen::VectorXi& patch_edge_is_init, Eigen::VectorXi& patch_vertex_is_init, int face, Eigen::MatrixXd& patch_vertices, Eigen::VectorXi& sharp_edge, Eigen::MatrixXi& tmp_sharp_edge, Eigen::VectorXi &edge_boundary_markers, Eigen::MatrixXi &tmp_edge_boundary_markers, Eigen::MatrixXd& V_orig, Eigen::VectorXi & vertex_is_fixed_orig, Eigen::VectorXi& new_mapped_indices_orig, Eigen::VectorXi& mesh_to_patch_indices) {
	if (sharp_edge(edge)) {

	}
	if (patch_edge_is_init[edge]) {
		return;
	}

	int start = EV(edge, 0);
	int end = EV(edge, 1);
	get_patch_vertex(start, face, patch_vertices, patch_vertex_is_init, V_orig, vertex_is_fixed_orig, new_mapped_indices_orig, mesh_to_patch_indices);
	get_patch_vertex(end, face, patch_vertices, patch_vertex_is_init, V_orig, vertex_is_fixed_orig, new_mapped_indices_orig, mesh_to_patch_indices);
	if (sharp_edge[edge]) {
		tmp_sharp_edge.conservativeResize(tmp_sharp_edge.rows() + 1, Eigen::NoChange);
		tmp_sharp_edge.bottomRows(1) << start, end;
	}
	tmp_edge_boundary_markers.conservativeResize(tmp_edge_boundary_markers.rows() + 1, Eigen::NoChange);
	tmp_edge_boundary_markers.bottomRows(1) << start, end, edge_boundary_markers[edge];
	patch_edge_is_init[edge] = 1;
}


void Patch::get_patch_vertex(int v_idx, int face, Eigen::MatrixXd& patch_vertices, Eigen::VectorXi& patch_vertex_is_init, Eigen::MatrixXd& V_orig, Eigen::VectorXi& vertex_is_fixed_orig, Eigen::VectorXi& new_mapped_indices_orig, Eigen::VectorXi& mesh_to_patch_indices) {
	if (patch_vertex_is_init[v_idx]) {
		return;
	}

	patch_vertices.conservativeResize(patch_vertices.rows() + 1, Eigen::NoChange);
	patch_vertices.bottomRows(1) = V_orig.row(v_idx);
	parent_vertices.conservativeResize(parent_vertices.rows() + 1);
	parent_vertices.tail(1) << v_idx;
	mesh_to_patch_indices[v_idx] = parent_vertices.rows() - 1;
	new_vertex_is_fixed.conservativeResize(new_vertex_is_fixed.rows() + 1);
	new_vertex_is_fixed.tail(1) << vertex_is_fixed_orig.row(v_idx);
	new_mapped_indices.conservativeResize(new_mapped_indices.rows() + 1);
	new_mapped_indices.tail(1) << -1;
	patch_vertex_is_init[v_idx] = 1;
}

void Patch::update_parent_vertex_positions(Eigen::MatrixXd& base_V) {
	for (int i = 0; i < mesh.V.rows(); i++) {
		base_V.row(parent_vertices[i]) = mesh.V.row(i);
	}
}

void Patch::update_patch_vertex_positions(Eigen::MatrixXd& base_V) {
	for (int i = 0; i < base_V.rows(); i++) {
		if (mesh_to_patch_indices[i] != -1) {
			mesh.V.row(mesh_to_patch_indices[i]) = base_V.row(i);
		}
	}
}
