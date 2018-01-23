#ifndef _Mesh_Extrusion_H_
#define _Mesh_Extrusion_H_
#include <Eigen/Core>
#include "Stroke.h"
#include "Mesh.h"
#include "SurfacePath.h"

class MeshExtrusion {

public:
	static void extrude_prepare(Stroke & base, SurfacePath & surface_path);
	static void extrude_main(Eigen::MatrixXd & V, Eigen::MatrixXi & F, Eigen::VectorXi & vertex_boundary_markers, Eigen::VectorXi & part_of_original_stroke, Eigen::VectorXi &new_mapped_indices, Eigen::VectorXi &sharp_edge, SurfacePath & surface_path, Stroke & stroke, Stroke & base, Eigen::Matrix4f model, Eigen::Matrix4f view, Eigen::Matrix4f proj, Eigen::Vector4f viewport);

private:
	static int prev_vertex_count;
	static int ID;

	static void generate_mesh(Mesh & m, Eigen::MatrixXd front_loop3D, Eigen::Vector3d center, Eigen::Vector3d x_vec, Eigen::Vector3d y_vec, Eigen::Vector3d offset, int nr_silhouette_vert, std::vector<int> loop_base_original_indices, std::vector<int> sil_original_indices);
	
	static void post_extrude_prepare_update_points(Stroke& stroke, SurfacePath& surface_path);
	static void post_extrude_main_update_points(Stroke & stroke, Eigen::MatrixXd new_positions);
	static void post_extrude_main_update_bindings(Stroke & base, SurfacePath & surface_path);

};

#endif