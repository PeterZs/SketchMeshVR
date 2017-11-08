#ifndef _Surface_Smoothing_H_
#define _Surface_Smoothing_H_
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <unordered_map>
#include "Mesh.h"


class SurfaceSmoothing {
public:
	static int prev_vertex_count;
	static Eigen::MatrixX3d vertex_normals;
	static Eigen::VectorXd curvatures;
	static std::unordered_map<int, Eigen::MatrixXd> precomputed_L; //Maps to the mesh ID
	static std::unordered_map<int, Eigen::MatrixX3d> precomputed_target_normals;
	static std::unordered_map<int, Eigen::SparseMatrix<double>> precompute_matrix_for_curvatures;
	static std::unordered_map<int, Eigen::SparseMatrix<double>> precompute_matrix_for_positions;


	static void smooth(Eigen::MatrixXd &V, Eigen::MatrixXi &F, Eigen::VectorXi &vertex_boundary_markers);
	static void clear_precomputed_matrices();
	static void smooth_main(Mesh &m);
	static Eigen::MatrixXd get_precomputed_L(Mesh &m);
	static void set_precomputed_L(Mesh & m, Eigen::MatrixXd & L);
	static Eigen::MatrixX3d get_precomputed_target_normals(Mesh & m);
	static void set_precomputed_target_normals(Mesh & m, Eigen::MatrixX3d target_normals);
	static Eigen::SparseMatrix<double> get_precompute_matrix_for_curvatures(Mesh & m);
	static void set_precompute_matrix_for_curvatures(Mesh & m, Eigen::SparseMatrix<double> curv);
	static Eigen::SparseMatrix<double> get_precompute_matrix_for_positions(Mesh & m);
	static void set_precompute_matrix_for_positions(Mesh & m, Eigen::SparseMatrix<double> pos);
	static Eigen::MatrixXd compute_laplacian_matrix(Mesh & m);
	static Eigen::MatrixX3d compute_target_normals(Mesh & m);
	static Eigen::VectorXd compute_target_curvatures(Mesh & m, Eigen::MatrixXd & L);
	static Eigen::MatrixX3d compute_vertex_laplacians(Mesh & m);
	static Eigen::SparseMatrix<double> compute_matrix_for_curvatures(Eigen::MatrixXd & L, int n, int n_fixed, Eigen::VectorXi vertex_boundary_markers, Eigen::VectorXd & laplacian_weights, Eigen::VectorXd & constraint_weights);
	static Eigen::SparseMatrix<double> compute_matrix_for_positions(Eigen::MatrixXd & L, Eigen::MatrixX3d laplacians, Mesh & m, Eigen::VectorXd laplacian_weights, Eigen::VectorXd constraint_weights);
	static Eigen::MatrixX3d get_target_laplacians(Mesh & m, Eigen::VectorXd target_curvatures, Eigen::MatrixX3d target_normals);
	static void move_vertices_to_satisfy_laplacians(Mesh & m, Eigen::MatrixXd & L, Eigen::MatrixX3d laplacians);
	static Eigen::VectorXd get_curvatures(Mesh & m);

};


#endif