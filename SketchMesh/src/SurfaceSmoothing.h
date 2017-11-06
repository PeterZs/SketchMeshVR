#ifndef _SurfaceSmoothing_
#define _SurfaceSmoothing
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <unordered_map>
#include "Mesh.h"

class SurfaceSmoothing {
public:
	static int prev_vertex_count;
	static unordered_map<Mesh, Eigen::MatrixXd> precomputed_L;
	static unordered_map<Mesh, Eigen::MatrixX3d> precomputed_target_normals;
	static unordered_map<Mesh, Eigen::MatrixXd> precompute_matrix_for_curvatures;
	static unordered_map<Mesh, Eigen::MatrixXd> precompute_matrix_for_positions;


	static void smooth(Eigen::MatrixXd &V, Eigen::MatrixXi &F, Eigen::VectorXi &vertex_boundary_markers);
	static void smooth_main(Mesh &m);
	static void clear_precomputed_matrices();
	static Eigen::MatrixXd get_precomputed_L(Mesh &m);
	static void set_precomputed_L(Mesh & m, Eigen::MatrixXd & L);
	static Eigen::MatrixX3d get_precomputed_target_normals(Mesh & m);
	static void set_precomputed_target_normals(Mesh & m, Eigen::MatrixX3d target_normals);
	static Eigen::MatrixXd compute_laplacian_matrix(Mesh & m);
	static Eigen::MatrixX3d compute_target_normals(Mesh & m);
	static Eigen::VectorXd compute_target_curvatures(Mesh & m, Eigen::MatrixXd & L);
};


#endif