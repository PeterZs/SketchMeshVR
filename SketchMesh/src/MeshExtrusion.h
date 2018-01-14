#ifndef _Mesh_Extrusion_H_
#define _Mesh_Extrusion_H_
#include <Eigen/Core>
#include "Stroke.h"
#include "Mesh.h"
#include "SurfacePath.h"

class MeshExtrusion {

public:
	static void extrude_prepare(Stroke & base, SurfacePath & surface_path);
	static void extrude_main(Eigen::MatrixXd &V, Eigen::MatrixXi &F, Eigen::VectorXi &vertex_boundary_markers, Eigen::VectorXi &part_of_original_stroke, SurfacePath & surface_path, Stroke & stroke);

private:
	static int prev_vertex_count;
	static int ID;

	static void generate_mesh(Mesh & m, Eigen::MatrixXd front_loop3D, Eigen::Vector3d center, Eigen::Vector3d x_vec, Eigen::Vector3d y_vec, Eigen::Vector3d offset, int nr_silhouette_vert, std::vector<int> loop_base_original_indices, std::vector<int> sil_original_indices);

};

#endif