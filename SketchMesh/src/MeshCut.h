#ifndef _Mesh_Cut_H_
#define _Mesh_Cut_H_
#include <Eigen/Core>
#include "Stroke.h"
#include "Mesh.h"
#include "SurfacePath.h"

class MeshCut {

public:
	static void cut(Eigen::MatrixXd & V, Eigen::MatrixXi & F, Eigen::VectorXi & vertex_boundary_markers, Eigen::VectorXi & part_of_original_stroke, Eigen::VectorXi & new_mapped_indices, Eigen::VectorXi & sharp_edge, Stroke & stroke);


private:
	static int prev_vertex_count;
	static int ID;
	static Eigen::MatrixXi EV, FE, EF;

	static void cut_main(Mesh& m, SurfacePath& surface_path, Stroke& stroke);
	static Eigen::MatrixXd resample_by_length_with_fixes(std::vector<int> path_vertices, Mesh & m, double unit_length);
	static void mesh_open_hole(Eigen::VectorXi & boundary_vertices, Mesh & m);


	static void post_cut_update_points(Stroke& stroke, SurfacePath& surface_path);

};


#endif

