#include <iostream>
#include <algorithm>
#include <unordered_map>
#include <igl/adjacency_matrix.h>
#include <igl/sum.h>
#include <igl/diag.h>
#include <igl/per_vertex_normals.h>
#include <igl/adjacency_list.h>
#include <numeric>
#include "SurfaceSmoothing.h"


using namespace std;
using namespace igl;

int SurfaceSmoothing::prev_vertex_count = -1;
Eigen::MatrixX3d SurfaceSmoothing::vertex_normals(0,3);
Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
Eigen::VectorXd SurfaceSmoothing::curvatures(0);
static Eigen::VectorXd xm, ym, zm;
int ID = 0;

//Map to the mesh IDs
std::unordered_map<int, Eigen::MatrixXd> SurfaceSmoothing::precomputed_L; 
std::unordered_map<int, Eigen::MatrixX3d> SurfaceSmoothing::precomputed_target_normals;
std::unordered_map<int, Eigen::SparseMatrix<double>> SurfaceSmoothing::precompute_matrix_for_curvatures;
std::unordered_map<int, Eigen::SparseMatrix<double>> SurfaceSmoothing::precompute_matrix_for_positions;

void SurfaceSmoothing::smooth(Eigen::MatrixXd &V, Eigen::MatrixXi &F, Eigen::VectorXi &vertex_boundary_markers) {
	if(prev_vertex_count != V.rows()) { //The mesh topology has changed, so we need to reset the precomputed matrices
		clear_precomputed_matrices();
		prev_vertex_count = V.rows();
	}
	Mesh m(V, F, vertex_boundary_markers, ID);
	ID++;

	smooth_main(m);
	V = m.V;
	F = m.F;
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

	igl::per_vertex_normals(m.V, m.F, vertex_normals);

	Eigen::MatrixX3d target_normals = get_precomputed_target_normals(m);
	if(target_normals.rows() == 0) {
		target_normals = compute_target_normals(m);
		set_precomputed_target_normals(m, target_normals);
	}
	Eigen::VectorXd target_curvatures = compute_target_curvatures(m, L);
	Eigen::MatrixX3d target_laplacians = get_target_laplacians(m, target_curvatures, target_normals);
	move_vertices_to_satisfy_laplacians(m, L, target_laplacians);

}

Eigen::MatrixXd SurfaceSmoothing::get_precomputed_L(Mesh &m) {
	return precomputed_L[m.ID]; //Will insert a new pair if there is none with the key m yet
}

void SurfaceSmoothing::set_precomputed_L(Mesh &m, Eigen::MatrixXd &L) {
	precomputed_L[m.ID] = L;
}

Eigen::MatrixX3d SurfaceSmoothing::get_precomputed_target_normals(Mesh &m) {
	return precomputed_target_normals[m.ID]; //Will insert a new pair if there is non with the key m yet
}

void SurfaceSmoothing::set_precomputed_target_normals(Mesh &m, Eigen::MatrixX3d target_normals) {
	precomputed_target_normals[m.ID] = target_normals;
}

Eigen::SparseMatrix<double> SurfaceSmoothing::get_precompute_matrix_for_curvatures(Mesh &m) {
	return precompute_matrix_for_curvatures[m.ID];
}

void SurfaceSmoothing::set_precompute_matrix_for_curvatures(Mesh &m, Eigen::SparseMatrix<double> curv) {
	precompute_matrix_for_curvatures[m.ID] = curv;
}

Eigen::SparseMatrix<double> SurfaceSmoothing::get_precompute_matrix_for_positions(Mesh &m) {
	return precompute_matrix_for_positions[m.ID];
}

void SurfaceSmoothing::set_precompute_matrix_for_positions(Mesh &m, Eigen::SparseMatrix<double> pos) {
	precompute_matrix_for_positions[m.ID] = pos;
}

//NOTE: TODO: might need to change the sign of U (due to different standards)
Eigen::MatrixXd SurfaceSmoothing::compute_laplacian_matrix(Mesh &m) {
	Eigen::SparseMatrix<double> A;
	adjacency_matrix(m.F, A);
	// sum each row 
	Eigen::SparseVector<double> Asum;
	sum(A, 1, Asum);
	// Convert row sums into diagonal of sparse matrix
	Eigen::SparseMatrix<double> Adiag;
	diag(Asum, Adiag);
	// Build uniform laplacian
	Eigen::SparseMatrix<double> U;
	U = A - Adiag;
	return (Eigen::MatrixXd) -U;
}

Eigen::MatrixX3d SurfaceSmoothing::compute_target_normals(Mesh &m) {
	Eigen::MatrixX3d normals;
	per_vertex_normals(m.V, m.F, PER_VERTEX_NORMALS_WEIGHTING_TYPE_UNIFORM, normals);
	//normals = vertex_normals;
	vector<vector<int>> neighbors;
	adjacency_list(m.F, neighbors);


	vector<vector<Eigen::Vector3d>> vs;
	for(int i = 0; i < m.V.rows(); i++) {
		vs.push_back(vector<Eigen::Vector3d>());
		for(int j = 0; j < neighbors[i].size(); j++) {
			vs[i].push_back(normals.row(neighbors[i][j]));
		}
	}


	int k = 5 * (int)sqrt(m.V.rows());

	for(int i = 0; i < k; i++) {
		double max_var = 0.0;

		for(int j = 0; j < m.V.rows(); j++) { //TODO: check if this is possible to implement without loop, and with only igl's normal methods
			if(m.vertex_boundary_markers[i]) {
				continue; //Don't change anything on boundary vertices
			}
			Eigen::Vector3d vec = Eigen::Vector3d::Zero();
			for(int l = 0; l < neighbors[j].size(); l++) {
				vec += vs[j][l];// vertex_normals.row(neighbors[j][l]);
			}
			vec.normalize();
			max_var = max(max_var, (vec - normals.row(j).transpose()).norm());
			normals.row(i) = vec;
		}

		if(max_var < 0.00001) {
			break;
		}
	}
	return normals; //TODO: NOTE THAT NORMALS AREN'T SET YET, BECAUSE WE DON'T HAVE ACCESS TO A VIEWER.DATA HERE
}

Eigen::VectorXd SurfaceSmoothing::compute_target_curvatures(Mesh &m, Eigen::MatrixXd &L) {
	double CONSTRAINT_WEIGHT = 0.1; //Weight for the constrained vertices (all), smaller = smoother
	Eigen::VectorXd laplacian_weights = Eigen::VectorXd::Constant(m.V.rows(), 1.0);
	Eigen::VectorXd constraint_weights = Eigen::VectorXd::Constant(m.V.rows(), CONSTRAINT_WEIGHT);
	Eigen::VectorXd current_curvature = get_curvatures(m);

	for(int i = 0; i < m.V.rows(); i++) {
		if(m.vertex_boundary_markers[i]) { //If the vertex is a boundary vertex
			constraint_weights[i] = CONSTRAINT_WEIGHT * 0.00001;
		}
	}

	int nr_boundary_vertices = (m.vertex_boundary_markers.array() > 0).count();
	Eigen::VectorXd b = Eigen::VectorXd::Zero(m.V.rows() + nr_boundary_vertices);
	int count = 0;
	for(int i = 0; i < m.V.rows(); i++) {
		if(m.vertex_boundary_markers[i] == 1) {
			//cout << current_curvature(i)*constraint_weights(i)*0.98 << endl;

			b[m.V.rows() + count] = max(0.0, current_curvature(i)*constraint_weights(i)*0.98);
			count++;
		}
	}

	Eigen::SparseMatrix<double> A = get_precompute_matrix_for_curvatures(m);
	Eigen::SparseMatrix<double> AT = A.transpose();
	if(A.rows() == 0 && A.cols() == 0) {
		A = compute_matrix_for_curvatures(L, m.V.rows(), nr_boundary_vertices, m.vertex_boundary_markers, laplacian_weights, constraint_weights);
		AT = A.transpose();
		solver.compute(AT*A);
		set_precompute_matrix_for_curvatures(m, A);
	}

	if((curvatures.size() == 0) || curvatures.size() != m.V.rows()){
		curvatures = Eigen::VectorXd(m.V.rows() + nr_boundary_vertices); //TODO: might want to change this line if we constrain curvature for all vertices (also non-boundary)
	}
	curvatures = solver.solve(AT*b);
	cout << AT*b << endl;

	curvatures.cwiseMax(0);
	return curvatures;
}

Eigen::MatrixX3d SurfaceSmoothing::compute_vertex_laplacians(Mesh &m) {
	vector<vector<int>> neighbors;
	adjacency_list(m.F, neighbors);
	Eigen::MatrixX3d laplacians(m.V.rows(),3);
	Eigen::Vector3d vec = Eigen::Vector3d::Zero();
	for(int i = 0; i < m.V.rows(); i++) {
		vec = Eigen::Vector3d::Zero();
		int nr_neighbors = neighbors[i].size();
		for(int j = 0; j < nr_neighbors; j++) {
			vec += m.V.row(neighbors[i][j]);
		}
		laplacians.row(i) = m.V.row(i) - (vec * (1.0 / nr_neighbors)).transpose();
	}
	return laplacians;
}

Eigen::VectorXd SurfaceSmoothing::get_curvatures(Mesh &m) {
	Eigen::VectorXd curvatures(m.V.rows());
	Eigen::MatrixX3d laplacians = compute_vertex_laplacians(m);
	//cout << "vertex normals " << laplacians << endl;
	for(int i = 0; i < m.V.rows(); i++) {
		curvatures[i] = laplacians.row(i).dot(vertex_normals.row(i).transpose()); //TODO: check that this is the correct way (not the standar definition of curvature, but what's used in java teddy)
	}
	return curvatures;
}

Eigen::SparseMatrix<double> SurfaceSmoothing::compute_matrix_for_curvatures(Eigen::MatrixXd &L, int n, int n_fixed, Eigen::VectorXi vertex_boundary_markers, Eigen::VectorXd &laplacian_weights, Eigen::VectorXd &constraint_weights) {
	Eigen::SparseMatrix<double> A(n + n_fixed, n);
	int count = 0;
	for(int i = 0; i < n; i++) {
		for(int j = 0; j < n; j++) {
			if(abs(L(i,j)* laplacian_weights[i]) > 0.000001) {
				A.insert(i, j) = L(i, j)*laplacian_weights[i];
			}
		}
		if(vertex_boundary_markers[i] == 1) { //Add extra constraints for boundary vertices
			A.insert(n + count, i) = constraint_weights[i]; //TODO: RIGHT NOW THIS DONE ONLY FOR THE ORIGINAL BOUNDARY VERTICES. WHEREAS IN TEDDY THEY DO IT FOR ALL VERTICES (SINCE CONSTRAIN_CURVARTURES_AT_BOUNDARY_ONLY IS FALSE)
			count++;
		}
	}
	return A;
}

Eigen::SparseMatrix<double> SurfaceSmoothing::compute_matrix_for_positions(Eigen::MatrixXd &L, Eigen::MatrixX3d laplacians, Mesh &m, Eigen::VectorXd laplacian_weights, Eigen::VectorXd constraint_weights) {
	int nr_boundary_vertices = (m.vertex_boundary_markers.array() > 0).count();
	Eigen::SparseMatrix<double> A(laplacians.rows() + nr_boundary_vertices, m.V.rows());

	int count = 0;
	for(int i = 0; i < laplacians.rows(); i++) {
		for(int j = 0; j < m.V.rows(); j++) {
			if(abs(L(i, j)*laplacian_weights(i)) > 0.000001) {
				A.insert(i, j) = L(i, j)*laplacian_weights[i];
			}
		}

		if(m.vertex_boundary_markers[i] == 1) {
			A.insert(laplacians.rows() + count, i) = constraint_weights[i];
			count++;
		}
	}
	return A;
}

Eigen::MatrixX3d SurfaceSmoothing::get_target_laplacians(Mesh &m, Eigen::VectorXd target_curvatures, Eigen::MatrixX3d target_normals) {
	Eigen::MatrixX3d target_laplacians(m.V.rows(), 3);
	for(int i = 0; i < m.V.rows(); i++) {
		target_laplacians.row(i) = vertex_normals.row(i) * target_curvatures[i]; //TODO: CHECK WHETHER TO USE VERTEX_NORMALS OR TARGET_NORMALS
	}
	return target_laplacians;
}

void SurfaceSmoothing::move_vertices_to_satisfy_laplacians(Mesh &m, Eigen::MatrixXd &L, Eigen::MatrixX3d laplacians) {
	double CONSTRAINT_WEIGHT = 10000;
	Eigen::VectorXd laplacian_weights = Eigen::VectorXd::Constant(m.V.rows(), 1, 1);
	Eigen::VectorXd constraint_weights = Eigen::VectorXd::Constant(m.V.rows(), 1, CONSTRAINT_WEIGHT);
	for(int i = 0; i < m.V.rows(); i++) {
		if(m.vertex_boundary_markers[i] == 1) {
			laplacian_weights[i] = 0.00001;
		}
	}

	Eigen::SparseMatrix<double> A = get_precompute_matrix_for_positions(m);
	Eigen::SparseMatrix<double> AT = A.transpose();
	if(A.rows() == 0 && A.cols() == 0) {
		A = compute_matrix_for_positions(L, laplacians, m, laplacian_weights, constraint_weights);
		AT = A.transpose();
		solver.compute(AT*A);
		set_precompute_matrix_for_positions(m, A);
	} 

	int nr_boundary_vertices = (m.vertex_boundary_markers.array() > 0).count();

	Eigen::VectorXd Bx(laplacians.rows() + nr_boundary_vertices); //laplacians has same number of entries as vertices in the mesh
	Eigen::VectorXd By(laplacians.rows() + nr_boundary_vertices);
	Eigen::VectorXd Bz(laplacians.rows() + nr_boundary_vertices);

	int count = 0;
	for(int i = 0; i < laplacians.rows(); i++) {
		Bx[i] = laplacians(i, 0) * laplacian_weights[i];
		By[i] = laplacians(i, 1) * laplacian_weights[i];
		Bz[i] = laplacians(i, 2) * laplacian_weights[i];

		if(m.vertex_boundary_markers[i] == 1) {
			Bx[laplacians.rows() + count] = m.V(i, 0) *constraint_weights[i];
			By[laplacians.rows() + count] = m.V(i, 1) *constraint_weights[i];
			Bz[laplacians.rows() + count] = m.V(i, 2) *constraint_weights[i];
			count++;
		}
	}

	//TODO: FROM HERE IT'S CHECKED AND CORRECT
	Eigen::VectorXd Vnewx;
	Eigen::VectorXd Vnewy;
	Eigen::VectorXd Vnewz;

	Vnewx = solver.solve(AT*Bx);
	Vnewy = solver.solve(AT*By);
	Vnewz = solver.solve(AT*Bz);

	for(int i = 0; i < m.V.rows(); i++) {
		if(m.vertex_boundary_markers[i] == 0) {
			//cout << Vnewx[i] << " " << Vnewy[i] << " " << Vnewz[i] << endl;
			m.V.row(i) << Vnewx[i], Vnewy[i], Vnewz[i];
		}
	}
}