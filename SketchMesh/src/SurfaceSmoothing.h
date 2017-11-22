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
    static double vertex_weight;
    static double edge_weight;
	static std::unordered_map<int, Eigen::MatrixXd> precomputed_L; //Maps to the mesh ID
	static std::unordered_map<int, Eigen::MatrixX3d> precomputed_target_normals;
	static std::unordered_map<int, Eigen::SparseMatrix<double>> precompute_matrix_for_LM_and_edges;
    static std::unordered_map<int, Eigen::SparseMatrix<double>> AT_for_LM_and_edges;
	static std::unordered_map<int, Eigen::SparseMatrix<double>> precompute_matrix_for_positions;
    static std::unordered_map<int, Eigen::SparseMatrix<double>> AT_for_positions;



	static void smooth(Eigen::MatrixXd &V, Eigen::MatrixXi &F, Eigen::VectorXi &vertex_boundary_markers);
	static void clear_precomputed_matrices();
	static void smooth_main(Mesh &m);
	static Eigen::MatrixXd get_precomputed_L(Mesh &m);
	static void set_precomputed_L(Mesh & m, Eigen::MatrixXd & L);
	static Eigen::SparseMatrix<double> get_precompute_matrix_for_LM_and_edges(Mesh & m);
	static void set_precompute_matrix_for_LM_and_edges(Mesh & m, Eigen::SparseMatrix<double> curv);
    static Eigen::SparseMatrix<double> get_AT_for_LM_and_edges(Mesh & m);
    static void set_AT_for_LM_and_edges(Mesh & m, Eigen::SparseMatrix<double> AT);

	static Eigen::SparseMatrix<double> get_precompute_matrix_for_positions(Mesh & m);
	static void set_precompute_matrix_for_positions(Mesh & m, Eigen::SparseMatrix<double> pos);
    static Eigen::SparseMatrix<double> get_AT_for_positions(Mesh & m);
    static void set_AT_for_positions(Mesh & m, Eigen::SparseMatrix<double> AT);
	static Eigen::MatrixXd compute_laplacian_matrix(Mesh & m);
	static Eigen::VectorXd compute_initial_curvature(Mesh & m);
	static Eigen::VectorXd compute_target_LMs(Mesh & m, Eigen::MatrixXd & L);
	static Eigen::VectorXd compute_target_edge_lengths(Mesh & m, Eigen::MatrixXd & L);
	static void compute_target_vertices(Mesh & m, Eigen::MatrixXd & L, Eigen::VectorXd & target_LMs, Eigen::VectorXd & target_edge_lengths);
	static Eigen::MatrixX3d compute_vertex_laplacians(Mesh & m);
	static Eigen::VectorXd get_curvatures(Mesh & m);

};


#endif
