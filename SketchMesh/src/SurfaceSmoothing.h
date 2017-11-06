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
	static unordered_map<Mesh, Eigen::MatrixXd> precomputed_target_normals;
	static unordered_map<Mesh, Eigen::MatrixXd> precompute_matrix_for_curvatures;
	static unordered_map<Mesh, Eigen::MatrixXd> precompute_matrix_for_positions;


	static void smooth(Eigen::MatrixXd &V, Eigen::MatrixXi &F);
	static void smooth_main(Mesh &m);
	static void clear_precomputed_matrices();
	static Eigen::MatrixXd get_precomputed_L(Mesh &m);
	static void set_precomputed_L(Mesh & m, Eigen::MatrixXd & L);
	static Eigen::MatrixXd compute_laplacian_matrix(Mesh & m);
};


#endif