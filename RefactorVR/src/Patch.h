#ifndef _PATCH_H_
#define _PATCH_H_

#include "Mesh.h"


class Patch {

public:
	Patch() : mesh(patch_vertices, patch_faces, new_vertex_boundary_markers, new_edge_boundary_markers, new_vertex_is_fixed, new_part_of_original_stroke, new_mapped_indices, new_sharp_edge, -1) {
		parent_vertices = Eigen::VectorXi();
		parent_faces = Eigen::VectorXi();
		new_vertex_boundary_markers = Eigen::VectorXi();
		new_edge_boundary_markers = Eigen::VectorXi();
		new_vertex_is_fixed = Eigen::VectorXi();
		new_part_of_original_stroke = Eigen::VectorXi();
		new_mapped_indices = Eigen::VectorXi();
		patch_vertex_is_init = Eigen::VectorXi();
		patch_edge_indices = Eigen::VectorXi();
		patch_edge_is_init = Eigen::VectorXi();
		new_sharp_edge = Eigen::VectorXi();
		patch_vertices = Eigen::MatrixXd();
		patch_faces = Eigen::MatrixXi();
		mesh_to_patch_indices = Eigen::VectorXi();
	;
	};

	static std::vector<Patch*> init_patches(Mesh& h);
	static void propagate_patch(Patch* patch, int face, Eigen::VectorXi& faces, std::vector<Patch*> &face_patch_map, Eigen::VectorXi& sharp_edge);
	void update_parent_vertex_positions(Eigen::MatrixXd& base_V);
	void update_patch_vertex_positions(Eigen::MatrixXd& base_V);
	//void update_patch_boundary_markers(Eigen::VectorXi& base_boundary_markers);
	Mesh get_mesh() { return mesh; };
	Mesh mesh;

private:
	void create_mesh_structure(Mesh& m, Eigen::VectorXi& faces);
	void get_patch_edge(int edge, Eigen::VectorXi & patch_edge_is_init, Eigen::VectorXi & patch_vertex_is_init, int face, Eigen::MatrixXd & patch_vertices, Eigen::VectorXi & sharp_edge, Eigen::VectorXi & new_sharp_edge, Eigen::VectorXi & edge_boundary_markers, Eigen::VectorXi & new_edge_boundary_markers, Eigen::MatrixXd & V_orig, Eigen::VectorXi & boundary_markers_orig, Eigen::VectorXi & vertex_is_fixed_orig, Eigen::VectorXi & part_of_original_orig, Eigen::VectorXi & new_mapped_indices_orig, Eigen::VectorXi & mesh_to_patch_indices);
	//void get_patch_edge(int edge, Eigen::VectorXi & patch_edge_is_init, Eigen::VectorXi & patch_vertex_is_init, int face, Eigen::MatrixXd & patch_vertices, Eigen::VectorXi & sharp_edge, Eigen::VectorXi & new_sharp_edge, Eigen::MatrixXd & V_orig, Eigen::VectorXi & boundary_markers_orig, Eigen::VectorXi & vertex_is_fixed_orig, Eigen::VectorXi & part_of_original_orig, Eigen::VectorXi & new_mapped_indices_orig, Eigen::VectorXi & mesh_to_patch_indices);
	//void get_patch_edge(int edge, Eigen::VectorXi & patch_edge_is_init, Eigen::VectorXi & patch_vertex_is_init, int face, Eigen::MatrixXd & patch_vertices, Eigen::VectorXi & sharp_edge, Eigen::VectorXi & new_sharp_edge, Eigen::MatrixXd & V_orig, Eigen::VectorXi & boundary_markers_orig, Eigen::VectorXi & part_of_original_orig, Eigen::VectorXi & new_mapped_indices_orig, Eigen::VectorXi & mesh_to_patch_indices);
//	void get_patch_vertex(int v_idx, int face, Eigen::MatrixXd & patch_vertices, Eigen::VectorXi & patch_vertex_is_init, Eigen::MatrixXd & V_orig, Eigen::VectorXi & boundary_markers_orig, Eigen::VectorXi & part_of_original_orig, Eigen::VectorXi & new_mapped_indices_orig, Eigen::VectorXi & mesh_to_patch_indices);
	void get_patch_vertex(int v_idx, int face, Eigen::MatrixXd & patch_vertices, Eigen::VectorXi & patch_vertex_is_init, Eigen::MatrixXd & V_orig, Eigen::VectorXi & boundary_markers_orig, Eigen::VectorXi & vertex_is_fixed_orig, Eigen::VectorXi & part_of_original_orig, Eigen::VectorXi & new_mapped_indices_orig, Eigen::VectorXi & mesh_to_patch_indices);

	Eigen::VectorXi parent_vertices, parent_faces, new_vertex_boundary_markers, new_edge_boundary_markers, new_vertex_is_fixed, new_part_of_original_stroke, new_mapped_indices;

	Eigen::MatrixXd patch_vertices; //Will form mesh.V
	Eigen::VectorXi patch_vertex_is_init;
	Eigen::MatrixXi patch_faces; //Will form mesh.f
	Eigen::VectorXi patch_edge_indices;
	Eigen::VectorXi patch_edge_is_init;
	Eigen::VectorXi new_sharp_edge; //Will form mesh.sharp_edge
	Eigen::VectorXi mesh_to_patch_indices;
};
#endif
