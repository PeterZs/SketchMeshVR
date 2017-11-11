#include <iostream>
#include <algorithm>
#include <unordered_map>
#include <igl/adjacency_matrix.h>
#include <igl/sum.h>
#include <igl/diag.h>
#include <igl/per_vertex_normals.h>
#include <igl/adjacency_list.h>
#include <igl/vertex_triangle_adjacency.h>
#include <numeric>
#include "SurfaceSmoothing.h"


using namespace std;
using namespace igl;

int SurfaceSmoothing::prev_vertex_count = -1;
Eigen::MatrixX3d SurfaceSmoothing::vertex_normals(0,3);
Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
Eigen::SparseLU<Eigen::SparseMatrix<double>> solver1;
Eigen::SparseLU<Eigen::SparseMatrix<double>> solver2;

Eigen::VectorXd SurfaceSmoothing::curvatures(0);
static Eigen::VectorXd xm, ym, zm;
int ID = 0;
bool firstDone = false;
int iteration = 0;

//Map to the mesh IDs
std::unordered_map<int, Eigen::MatrixXd> SurfaceSmoothing::precomputed_L; 
std::unordered_map<int, Eigen::MatrixX3d> SurfaceSmoothing::precomputed_target_normals;
std::unordered_map<int, Eigen::SparseMatrix<double>> SurfaceSmoothing::precompute_matrix_for_LM_and_edges;
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
	iteration++;
}

void SurfaceSmoothing::clear_precomputed_matrices() {
	cout << "Clearing precomputed matrices" << endl;
	precomputed_L.clear(); //NOTE: as opposed to java's Hashmap.clear, this reduces the size of the unordered_map to 0

	for(auto it = precompute_matrix_for_LM_and_edges.begin(); it != precompute_matrix_for_LM_and_edges.end(); ++it) {
		it->second.resize(0, 0); //Resize the referenced Eigen::MatrixXd to size 0
	}
	precompute_matrix_for_LM_and_edges.clear();

	for(auto it = precompute_matrix_for_positions.begin(); it != precompute_matrix_for_positions.end(); ++it) {
		it->second.resize(0, 0); //Resize the referenced Eigen::MatrixXd to size 0
	}
	precompute_matrix_for_positions.clear(); 

	precomputed_target_normals.clear();//TODO1: do we still need this?
}

void SurfaceSmoothing::smooth_main(Mesh &m) {
	Eigen::MatrixXd L = get_precomputed_L(m);
	if(L.rows()==0 && L.cols()==0) {
		L = compute_laplacian_matrix(m);
		set_precomputed_L(m, L);
	}

	igl::per_vertex_normals(m.V, m.F, vertex_normals);

	Eigen::VectorXd target_LMs = compute_target_LMs(m, L);
	Eigen::VectorXd target_edge_lengths = compute_target_edge_lengths(m, L);
	Eigen::MatrixXd target_vertices = compute_target_vertices(m, L, target_LMs, target_edge_lengths);


	//TODO1
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
} //TODO1

void SurfaceSmoothing::set_precomputed_L(Mesh &m, Eigen::MatrixXd &L) {
	precomputed_L[m.ID] = L;
} //TODO1

Eigen::MatrixX3d SurfaceSmoothing::get_precomputed_target_normals(Mesh &m) {
	return precomputed_target_normals[m.ID]; //Will insert a new pair if there is non with the key m yet
} //TODO1

void SurfaceSmoothing::set_precomputed_target_normals(Mesh &m, Eigen::MatrixX3d target_normals) {
	precomputed_target_normals[m.ID] = target_normals;
} //TODO1

Eigen::SparseMatrix<double> SurfaceSmoothing::get_precompute_matrix_for_LM_and_edges(Mesh &m) {
	return precompute_matrix_for_LM_and_edges[m.ID];
}

void SurfaceSmoothing::set_precompute_matrix_for_LM_and_edges(Mesh &m, Eigen::SparseMatrix<double> LM_edge) {
	precompute_matrix_for_LM_and_edges[m.ID] = LM_edge;
}

Eigen::SparseMatrix<double> SurfaceSmoothing::get_precompute_matrix_for_positions(Mesh &m) {
	return precompute_matrix_for_positions[m.ID];
}

void SurfaceSmoothing::set_precompute_matrix_for_positions(Mesh &m, Eigen::SparseMatrix<double> pos) {
	precompute_matrix_for_positions[m.ID] = pos;
}

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
//TODO1
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

Eigen::VectorXd SurfaceSmoothing::compute_target_LMs(Mesh &m, Eigen::MatrixXd &L) {
	Eigen::SparseMatrix<double> A = get_precompute_matrix_for_LM_and_edges(m);
	Eigen::SparseMatrix<double> AT = A.transpose();
	if(A.rows() == 0 && A.cols() == 0) { //We haven't set up A for this topology yet
		//A.block(0, 0, L.rows(), L.cols()) = L;
		for(int i = 0; i < m.V.rows(); i++) {
			for(int j = 0; j < m.V.rows(); j++) {
				A.insert(i, j) = L(i, j);
			}
			if(m.vertex_boundary_markers[i] == 1) { //Constrain only the boundary in the first iteration
				A.insert(m.V.rows() + i, i) = 1; //For target LM'/edge'
			}
		}
		AT = A.transpose();
		solver1.compute(AT*A);
		set_precompute_matrix_for_LM_and_edges(m, A);
	}
	if(iteration == 1) { //Constrain all vertices after the first iteration
		for(int i = 0; i < m.V.rows(); i++) {
			if(m.vertex_boundary_markers[i] == 0) {
				A.insert(m.V.rows() + i, i) = 1; 
			}
		}
		AT = A.transpose();
		solver1.compute(AT*A);
		set_precompute_matrix_for_LM_and_edges(m, A);
	}

	Eigen::VectorXd current_curvatures = get_curvatures(m); //TODO1: email about this being the right curvature (for all iterations) or not
	Eigen::VectorXd b = Eigen::VectorXd::Zero(m.V.rows()*2);
	for(int i = 0; i < m.V.rows(); i++) {
		if(iteration == 0) {
			if(m.vertex_boundary_markers[i] == 1) {
				b[m.V.rows() + i] = current_curvatures[i];
			}
		} else {
			b[m.V.rows() + i] = current_curvatures[i];
		}
	}
	Eigen::VectorXd target_LM;
	target_LM = solver1.solve(AT*b);

	return target_LM;
}

Eigen::VectorXd SurfaceSmoothing::compute_target_edge_lengths(Mesh &m, Eigen::MatrixXd &L) {
	Eigen::SparseMatrix<double> A = get_precompute_matrix_for_LM_and_edges(m);
	Eigen::SparseMatrix<double> AT = A.transpose();
	if(A.rows() == 0 && A.cols() == 0) { //We haven't set up A for this topology yet
		//THIS SHOULD NEVER HAPPEN, ASSUMING THE TARGET LMS ALWAYS GET COMPUTED BEFORE THE TARGET EDGE LENGTHS
		cout << "this shouldn't happen" << endl;
		//A.block(0, 0, L.rows(), L.cols()) = L;
		for(int i = 0; i < m.V.rows(); i++) {
			for(int j = 0; j < m.V.rows(); j++) {
				A.insert(i, j) = L(i, j);
			}
			if(m.vertex_boundary_markers[i] == 1) { //Constrain only the boundary in the first iteration
				A.insert(m.V.rows() + i, i) = 1; //For target LM'/edge'
			}
		}
		AT = A.transpose();
		solver1.compute(AT*A);
		set_precompute_matrix_for_LM_and_edges(m, A);
	}
	if(iteration == 1) { //Constrain all vertices after the first iteration
		cout << A << endl;
		cout << " A should already be filled everywhere" << endl;
		for(int i = 0; i < m.V.rows(); i++) {
			if(m.vertex_boundary_markers[i] == 0) {
				A.insert(m.V.rows() + i, i) = 1;
			}
		}
		AT = A.transpose();
		solver1.compute(AT*A);
		set_precompute_matrix_for_LM_and_edges(m, A);
	}

	vector<vector<int>> neighbors;
	adjacency_list(m.F, neighbors);
	Eigen::VectorXd b = Eigen::VectorXd::Zero(m.V.rows() * 2);
	for(int i = 0; i < m.V.rows(); i++) {
		if(iteration == 0) {
			if(m.vertex_boundary_markers[i] == 1) {
				double cum_edge_length = 0.0;
				for(int j = 0; j < neighbors[i].size(); j++) {
					cum_edge_length += (m.V.row(neighbors[i][j]) - m.V.row(i)).norm();
				}
				b[m.V.rows() + i] = (double) (cum_edge_length / neighbors[i].size());
			}
		} else {
			double cum_edge_length = 0.0;
			for(int j = 0; j < neighbors[i].size(); j++) {
				cum_edge_length += (m.V.row(neighbors[i][j]) - m.V.row(i)).norm();
			}
			b[m.V.rows() + i] = (double) (cum_edge_length /neighbors[i].size());
		}
	}

	Eigen::VectorXd target_edge_lengths;
	target_edge_lengths = solver1.solve(AT*b);

	return target_edge_lengths;
}

Eigen::MatrixX3d SurfaceSmoothing::compute_target_vertices(Mesh &m, Eigen::MatrixXd &L, Eigen::VectorXd &target_LMs, Eigen::VectorXd &target_edge_lengths) {
	double vertex_weight = 100, edge_weight = 0.001;
	vector<vector<int>> neighbors;
	adjacency_list(m.F, neighbors);
	int no_boundary_vertices = (m.vertex_boundary_markers.array() > 0).count();
	int no_boundary_adjacent_vertices = 0;

	for(int i = 0; i < m.V.rows(); i++) {
		if(m.vertex_boundary_markers[i] == 1) {
			no_boundary_adjacent_vertices += neighbors[i].size();
		}
	}

	Eigen::SparseMatrix<double> A = get_precompute_matrix_for_positions(m);
	Eigen::SparseMatrix<double> AT = A.transpose();
	if(A.rows() == 0 && A.cols() == 0) { //We haven't set up A for this topology yet
		//A.block(0, 0, L.rows(), L.cols()) = L;
		int count = 0, count2=0;
		for(int i = 0; i < m.V.rows(); i++) {
			for(int j = 0; j < m.V.rows(); j++) {
				A.insert(i, j) = L(i, j);
			}
			if(m.vertex_boundary_markers[i] == 1) { //Constrain only the boundary in the first iteration
				A.insert(m.V.rows() + count, i) = vertex_weight; //For target LM'/edge'
				count++;

				for(int j = 0; j < neighbors[i].size(); j++) {
					A.insert(m.V.rows() + no_boundary_vertices + count2, i) = edge_weight;
					A.insert(m.V.rows() + no_boundary_vertices + count2, neighbors[i][j]) = -edge_weight;
				}
			}
		}
		AT = A.transpose();
		solver2.compute(AT*A);
		set_precompute_matrix_for_LM_and_edges(m, A);
	}

	Eigen::VectorXd bx = Eigen::VectorXd::Zero(m.V.rows() + no_boundary_vertices + no_boundary_adjacent_vertices);
	Eigen::VectorXd by = Eigen::VectorXd::Zero(m.V.rows() + no_boundary_vertices + no_boundary_adjacent_vertices);
	Eigen::VectorXd bz = Eigen::VectorXd::Zero(m.V.rows() + no_boundary_vertices + no_boundary_adjacent_vertices);


	Eigen::VectorXd doubleAreas;
	igl::doublearea(m.V, m.F, doubleAreas);
	vector<vector<int>> VF, VFi;
	igl::vertex_triangle_adjacency(m.V.rows(), m.F, VF, VFi);
	int count = 0, count2 = 0;
	for(int i = 0; i < m.V.rows(); i++) {
		double vertex_area = 0.0;
		for(int j = 0; j < neighbors[i].size(); j++) {
			vertex_area += doubleAreas(VF[i][j]); //Get the area of the adjacent face
		}
		vertex_area /= 3.0;

		Eigen::Vector3d delta = vertex_area * target_LMs[i] * vertex_normals.row(i);
		bx[i] = delta(0);
		by[i] = delta(1);
		bz[i] = delta(2);

		if(m.vertex_boundary_markers[i] == 1) {
			for(int j = 0; j < neighbors[i].size(); j++) {
				Eigen::Vector3d vec = (m.V.row(i) - m.V.row(neighbors[i][j]));
				vec.normalize();
				Eigen::Vector3d edge_vec = (target_edge_lengths[i] - target_edge_lengths[neighbors[i][j]]) * vec;

				bx[m.V.rows() + no_boundary_vertices + count] = edge_vec(0) * edge_weight;
				by[m.V.rows() + no_boundary_vertices + count] = edge_vec(1) * edge_weight;
				bz[m.V.rows() + no_boundary_vertices + count] = edge_vec(2) * edge_weight;
				count++;

				bx[m.V.rows() + count2] = m.V(i, 0) * vertex_weight;
				by[m.V.rows() + count2] = m.V(i, 1) * vertex_weight;
				bz[m.V.rows() + count2] = m.V(i, 2) * vertex_weight;
				count2++;
			}
		}
	}

	Eigen::VectorXd Vnewx;
	Eigen::VectorXd Vnewy;
	Eigen::VectorXd Vnewz;

	Vnewx = solver2.solve(AT*bx);
	Vnewy = solver2.solve(AT*by);
	Vnewz = solver2.solve(AT*bz);

	//TODO: check if we should do this for all vertices instead
	for(int i = 0; i < m.V.rows(); i++) {
		if(m.vertex_boundary_markers[i] == 0) {
			//cout << Vnewx[i] << " " << Vnewy[i] << " " << Vnewz[i] << endl;
			m.V.row(i) << Vnewx[i], Vnewy[i], Vnewz[i];
		}
	}

	Eigen::MatrixX3d newVertices(m.V.rows(), 3);
	newVertices.col(0) << Vnewx;
	newVertices.col(1) << Vnewy;
	newVertices.col(2) << Vnewz;

	return newVertices;
}
//TODO1
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
		if(firstDone) { //TODO: CHECK IF WE SHOULD ALSO SET TARGET CURVATURES FOR NON-BOUNDARY VERTICES AFTER THE FIRST ROUND
			b[i] = current_curvature(i)*constraint_weights(i);
			if(m.vertex_boundary_markers[i] == 1) {
				b[m.V.rows() + count] = current_curvature(i)*constraint_weights(i)*0.98;
				count++;
			}
		} else {
			if(m.vertex_boundary_markers[i] == 1) {
				b[m.V.rows() + count] = current_curvature(i)*constraint_weights(i)*0.98;
				count++;
				firstDone = true;
			}
		}
	}

	Eigen::SparseMatrix<double> A = get_precompute_matrix_for_LM_and_edges(m);
	Eigen::SparseMatrix<double> AT = A.transpose();
	if(A.rows() == 0 && A.cols() == 0) {
		A = compute_matrix_for_curvatures(L, m.V.rows(), nr_boundary_vertices, m.vertex_boundary_markers, laplacian_weights, constraint_weights);
		AT = A.transpose();
		solver.compute(AT*A);
		set_precompute_matrix_for_LM_and_edges(m, A);
	}

	if((curvatures.size() == 0) || curvatures.size() != m.V.rows()){
		curvatures = Eigen::VectorXd(m.V.rows() + nr_boundary_vertices); //TODO: might want to change this line if we constrain curvature for all vertices (also non-boundary)
	}
	curvatures = solver.solve(AT*b);

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
	for(int i = 0; i < m.V.rows(); i++) {
		curvatures[i] = laplacians.row(i).dot(vertex_normals.row(i).transpose()); //TODO: check that this is the correct way (not the standar definition of curvature, but what's used in java teddy)
	}
	return curvatures;
}
//TODO1
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
//TODO1
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
		target_laplacians.row(i) = vertex_normals.row(i) * target_curvatures[i]; //TODO: CHECK WHETHER TO USE VERTEX_NORMALS OR TARGET_NORMALS. AND WHETHER TO INCLUDE THE AREA CORRESPONDING TO THE VERTEX
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