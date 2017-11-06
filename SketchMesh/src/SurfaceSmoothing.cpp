#include "SurfaceSmoothing.h"
#include <iostream>
#include <igl/adjacency_matrix.h>
#include <igl/sum.h>
#include <igl/diag.h>
using namespace std;
using namespace igl;

int SurfaceSmoothing::prev_vertex_count = -1;

void SurfaceSmoothing::smooth(Eigen::MatrixXd &V, Eigen::MatrixXi &F) {
	if(prev_vertex_count != V.rows()) { //The mesh topology has changed, so we need to reset the precomputed matrices
		clear_precomputed_matrices();
		prev_vertex_count = V.rows();
	}
	Mesh m(V, F);
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
}

Eigen::MatrixXd SurfaceSmoothing::get_precomputed_L(Mesh &m) {
	return precomputed_L[m]; //Will insert a new pair if there is none with the key m yet
}

void SurfaceSmoothing::set_precomputed_L(Mesh &m, Eigen::MatrixXd &L) {
	precomputed_L[m] = L;
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
