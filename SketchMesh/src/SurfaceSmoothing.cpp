#include "SurfaceSmoothing.h"
#include <iostream>
#include <algorithm>
#include <igl/adjacency_matrix.h>
#include <igl/sum.h>
#include <igl/diag.h>
#include <igl/per_vertex_normals.h>
#include <igl/adjacency_list.h>
using namespace std;
using namespace igl;

int SurfaceSmoothing::prev_vertex_count = -1;
Eigen::MatrixX3d vertex_normals;


void SurfaceSmoothing::smooth(Eigen::MatrixXd &V, Eigen::MatrixXi &F, Eigen::VectorXi &vertex_boundary_markers) {
	if(prev_vertex_count != V.rows()) { //The mesh topology has changed, so we need to reset the precomputed matrices
		clear_precomputed_matrices();
		prev_vertex_count = V.rows();
	}
	Mesh m(V, F, vertex_boundary_markers);
	smooth_main(m);
}

void SurfaceSmoothing::clear_precomputed_matrices() {
	cout << "Clearing precomputed matrices" << endl;
	precomputed_L.clear(); //NOTE: as opposed to java's Hashmap.clear, this reduces the size of the unordered_map to 0

	for(auto it = precompute_matrix_for_curvatures.begin(); it != precompute_matrix_for_curvatures.end(); ++it) {
		it->second.resize(0, 0); //Resize the referenced Eigen::MatrixXd to size 0
	}
	precompute_matrix_for_curvatures.clear();

	for(auto it = precompute_matrix_for_positions.begin(); it != precompute_matrix_for_positions.end(); ++it) {
		it->second.resize(0, 0); //Resize the referenced Eigen::MatrixXd to size 0
	}
	precompute_matrix_for_positions.clear();

	precomputed_target_normals.clear();
}

void SurfaceSmoothing::smooth_main(Mesh &m) {
	Eigen::MatrixXd L = get_precomputed_L(m);
	if(L.rows()==0 && L.cols()==0) {
		L = compute_laplacian_matrix(m);
		set_precomputed_L(m, L);
	}

	Eigen::MatrixX3d target_normals = get_precomputed_target_normals(m);
	if(target_normals.rows() == 0 && target_normals.cols() == 0) {
		target_normals = compute_target_normals(m);
		set_precomputed_target_normals(m, target_normals);
	}
	Eigen::VectorXd target_curvatures = compute_target_curvatures(m, L);
	Eigen::VectorXcd target_laplacians = get_target_laplacians(m, target_curvatures, target_normals);
	move_vertices_to_satisfy_laplacians(m, L, target_laplacians);

}

Eigen::MatrixXd SurfaceSmoothing::get_precomputed_L(Mesh &m) {
	return precomputed_L[m]; //Will insert a new pair if there is none with the key m yet
}

void SurfaceSmoothing::set_precomputed_L(Mesh &m, Eigen::MatrixXd &L) {
	precomputed_L[m] = L;
}

Eigen::MatrixX3d SurfaceSmoothing::get_precomputed_target_normals(Mesh &m) {
	return precomputed_target_normals[m]; //Will insert a new pair if there is non with the key m yet
}

void SurfaceSmoothing::set_precomputed_target_normals(Mesh &m, Eigen::MatrixX3d target_normals) {
	precomputed_target_normals[m] = target_normals;
}

//NOTE: TODO: might need to change the sign of U (due to different standards)
Eigen::MatrixXd SurfaceSmoothing::compute_laplacian_matrix(Mesh &m) {
	Eigen::SparseMatrix<double> A;
	adjacency_matrix(m.F, A);
	Eigen::SparseMatrix<double> A;
	// sum each row 
	Eigen::SparseVector<double> Asum;
	sum(A,1,Asum);
	// Convert row sums into diagonal of sparse matrix
	Eigen::SparseMatrix<double> Adiag;
	diag(Asum,Adiag);
	// Build uniform laplacian
	Eigen::SparseMatrix<double> U;
	U = A-Adiag;
	return (Eigen::MatrixXd) U;
}

Eigen::MatrixX3d SurfaceSmoothing::compute_target_normals(Mesh &m) {
	Eigen::MatrixX3d normals;
	per_vertex_normals(m.V, m.F, PER_VERTEX_NORMALS_WEIGHTING_TYPE_UNIFORM, vertex_normals);
	normals = vertex_normals;
	vector<vector<int>> neighbors;
	adjacency_list(m.F, neighbors);
	int k = 5 * (int)sqrt(m.V.rows());

	for(int i = 0; i < k; i++) {
		double max_var = 0;

		for(int j = 0; j < m.V.rows(); j++) { //TODO: check if this is possible to implement without loop, and with only igl's normal methods
			if(m.vertex_boundary_markers[i]) {
				continue; //Don't change anything on boundary vertices
			}
			Eigen::Vector3d vec = Eigen::Vector3d::Zero();
			for(int l = 0; l < neighbors[j].size(); l++) {
				vec += vertex_normals.row(neighbors[j][l]);
			}
			vec.normalize();
			max_var = max(max_var, (vec - normals.row(j)).norm());

			normals.row(i) = vec;
		}


		if(max_var < 0.00001) {
			break;
		}
	}
	return normals; //TODO: NOTE THAT NORMALS AREN'T SET YET, BECAUSE WE DON'T HAVE ACCESS TO A VIEWER.DATA HERE
}


Eigen::VectorXd compute_target_curvatures(Mesh &m, Eigen::MatrixXd &L) {
	double CONSTRAINT_WEIGHT = 0.1; //Weight for the constrained vertices (all), smaller = smoother
	Eigen::VectorXd laplacian_weights = Eigen::VectorXd::Constant(m.V.rows(), 1.0);
	Eigen::VectorXd constraint_weights = Eigen::VectorXd::Constant(m.V.rows(), CONSTRAINT_WEIGHT);
	for(int i = 0; i < m.V.rows(); i++) {
		if(m.vertex_boundary_markers[i]) { //If the vertex is a boundary vertex
			constraint_weights[i] = CONSTRAINT_WEIGHT * 0.00001;
		}
		current_curvature[i] = 
	}
}

double SurfaceSmoothing::get_curvature