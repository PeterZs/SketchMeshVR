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
Eigen::VectorXd initial_curvature;

//Map to the mesh IDs
std::unordered_map<int, Eigen::MatrixXd> SurfaceSmoothing::precomputed_L; 
std::unordered_map<int, Eigen::MatrixX3d> SurfaceSmoothing::precomputed_target_normals;
std::unordered_map<int, Eigen::SparseMatrix<double>> SurfaceSmoothing::precompute_matrix_for_LM_and_edges;
std::unordered_map<int, Eigen::SparseMatrix<double>> SurfaceSmoothing::precompute_matrix_for_positions;

void SurfaceSmoothing::smooth(Eigen::MatrixXd &V, Eigen::MatrixXi &F, Eigen::VectorXi &vertex_boundary_markers) {
	if(prev_vertex_count != V.rows()) { //The mesh topology has changed, so we need to reset the precomputed matrices
		clear_precomputed_matrices();
		prev_vertex_count = V.rows();
		iteration = 0;
		cout << " reset" << endl;
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
	igl::per_vertex_normals(m.V, m.F, vertex_normals); //TODO: might want to use PER_VERTEX_NORMALS_WEIGHTING_TYPE_UNIFORM
	Eigen::VectorXd target_LMs = compute_target_LMs(m, L);
	Eigen::VectorXd target_edge_lengths = compute_target_edge_lengths(m, L);
	Eigen::MatrixX3d target_vertices = compute_target_vertices(m, L, target_LMs, target_edge_lengths);
}

Eigen::MatrixXd SurfaceSmoothing::get_precomputed_L(Mesh &m) {
	return precomputed_L[m.ID]; //Will insert a new pair if there is none with the key m yet
} 

void SurfaceSmoothing::set_precomputed_L(Mesh &m, Eigen::MatrixXd &L) {
	precomputed_L[m.ID] = L;
} 

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

Eigen::VectorXd SurfaceSmoothing::compute_initial_curvature(Mesh &m) {
	int no_boundary_vertices = (m.vertex_boundary_markers.array() > 0).count();
	Eigen::VectorXd initial_curvatures(no_boundary_vertices);
	initial_curvatures[0] = (m.V.row(0) - 0.5 * (m.V.row(no_boundary_vertices - 1) + m.V.row((0 + 1) % no_boundary_vertices))).norm();

	for(int i = 1; i <no_boundary_vertices; i++) {
		initial_curvatures[i] = (m.V.row(i) - 0.5 * (m.V.row(i - 1) + m.V.row((i + 1) % no_boundary_vertices))).norm();
	}
	return initial_curvatures;
}

Eigen::VectorXd SurfaceSmoothing::compute_target_LMs(Mesh &m, Eigen::MatrixXd &L) {
	Eigen::SparseMatrix<double> A = get_precompute_matrix_for_LM_and_edges(m);
	Eigen::SparseMatrix<double> AT = A.transpose();
	if(A.rows() == 0 && A.cols() == 0) { //We haven't set up A for this topology yet
		A = Eigen::SparseMatrix<double>(m.V.rows()*2, m.V.rows());
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

	if(iteration == 0) {
		initial_curvature = compute_initial_curvature(m);
	}

	Eigen::VectorXd current_curvatures = get_curvatures(m); 
	Eigen::VectorXd b = Eigen::VectorXd::Zero(m.V.rows()*2);
	int count = 0;
	for(int i = 0; i < m.V.rows(); i++) {
		if(iteration == 0) {
			if(m.vertex_boundary_markers[i] == 1) {
				b[m.V.rows() + i] = initial_curvature[i];
				count++;
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
		A = Eigen::SparseMatrix<double>(m.V.rows() + no_boundary_vertices + no_boundary_adjacent_vertices, m.V.rows());
		int count = 0, count2 = 0;
		for(int i = 0; i < m.V.rows(); i++) {
			for(int j = 0; j < m.V.rows(); j++) {
				A.insert(i, j) = L(i, j);
			}
			if(m.vertex_boundary_markers[i] == 1) { //Constrain only the boundary LMs and edges
				A.insert(m.V.rows() + count, i) = vertex_weight; //For target LM'/edge'
				count++;

				for(int j = 0; j < neighbors[i].size(); j++) {
					A.insert(m.V.rows() + no_boundary_vertices + count2, i) = edge_weight;
					A.insert(m.V.rows() + no_boundary_vertices + count2, neighbors[i][j]) = -edge_weight;
					count2++;
				}
			}
		}
		AT = A.transpose();
		solver2.compute(AT*A);
		set_precompute_matrix_for_positions(m, A);
	}

	Eigen::VectorXd bx = Eigen::VectorXd::Zero(m.V.rows() + no_boundary_vertices + no_boundary_adjacent_vertices);
	Eigen::VectorXd by = Eigen::VectorXd::Zero(m.V.rows() + no_boundary_vertices + no_boundary_adjacent_vertices);
	Eigen::VectorXd bz = Eigen::VectorXd::Zero(m.V.rows() + no_boundary_vertices + no_boundary_adjacent_vertices);

	Eigen::VectorXd doubleAreas;
	igl::doublearea(m.V, m.F, doubleAreas);
	double totalArea = doubleAreas.sum() / 2.0;
	
	
	vector<vector<int>> VF, VFi;
	igl::vertex_triangle_adjacency(m.V.rows(), m.F, VF, VFi);
	Eigen::VectorXd vertex_areas(m.V.rows());
	for(int i = 0; i < m.V.rows(); i++) {
		double vertex_area = 0.0;
		for(int j = 0; j < VF[i].size(); j++) {
			vertex_area += doubleAreas(VF[i][j]); //Get the area of the adjacent face
		}
		vertex_areas[i] = vertex_area / 6.0; //Divide by 2 because of "double" area, and then by 3 because every vertex shares a face with 2 other vertices
	}
	double min_area = vertex_areas.minCoeff();

	int count = 0, count2 = 0;
	for(int i = 0; i < m.V.rows(); i++) {
		Eigen::Vector3d delta = vertex_areas[i]/min_area * target_LMs[i] * vertex_normals.row(i);
		bx[i] = delta(0);
		by[i] = delta(1);
		bz[i] = delta(2);

		if(m.vertex_boundary_markers[i] == 1) {
			bx[m.V.rows() + count2] = m.V(i, 0) * vertex_weight;
			by[m.V.rows() + count2] = m.V(i, 1) * vertex_weight;
			bz[m.V.rows() + count2] = m.V(i, 2) * vertex_weight;
			count2++;

			for(int j = 0; j < neighbors[i].size(); j++) {
				Eigen::Vector3d vec = (m.V.row(i) - m.V.row(neighbors[i][j]));
				vec.normalize();
				Eigen::Vector3d edge_vec = (target_edge_lengths[i] + target_edge_lengths[neighbors[i][j]])/2.0 * vec;

				bx[m.V.rows() + no_boundary_vertices + count] = edge_vec(0) * edge_weight;
				by[m.V.rows() + no_boundary_vertices + count] = edge_vec(1) * edge_weight;
				bz[m.V.rows() + no_boundary_vertices + count] = edge_vec(2) * edge_weight;
				count++;
			}
		}
	}

	Eigen::VectorXd Vnewx;
	Eigen::VectorXd Vnewy;
	Eigen::VectorXd Vnewz;

	Vnewx = solver2.solve(AT*bx);
	Vnewy = solver2.solve(AT*by);
	Vnewz = solver2.solve(AT*bz);

	for(int i = 0; i < m.V.rows(); i++) {
		if(m.vertex_boundary_markers[i] == 0) {
			m.V.row(i) << Vnewx[i], Vnewy[i], Vnewz[i];
		}
	}
	Eigen::MatrixX3d newVertices(m.V.rows(), 3);
	newVertices.col(0) << Vnewx;
	newVertices.col(1) << Vnewy;
	newVertices.col(2) << Vnewz;
	cout << m.V << endl << endl;
	return newVertices;
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