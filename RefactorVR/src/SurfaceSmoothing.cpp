#include <iostream>
#include <algorithm>
#include <unordered_map>
#include <igl/adjacency_matrix.h>
#include <igl/sum.h>
#include <igl/diag.h>
#include <igl/per_vertex_normals.h>
#include <igl/adjacency_list.h>
#include <igl/edge_topology.h>
#include "SurfaceSmoothing.h"
#include "Patch.h"

using namespace std;
using namespace igl;

typedef Eigen::Triplet<double> T;

int SurfaceSmoothing::prev_vertex_count = -1;
static int prev_mesh_ID = -1;
Eigen::MatrixX3d SurfaceSmoothing::vertex_normals(0, 3);
//Eigen::SparseLU<Eigen::SparseMatrix<double>> solver1;
//Eigen::SparseLU<Eigen::SparseMatrix<double>> solver2;
std::vector<Eigen::SparseLU<Eigen::SparseMatrix<double>>*> solver1_array;
std::vector<Eigen::SparseLU<Eigen::SparseMatrix<double>>*> solver2_array;

Eigen::VectorXd initial_curvature;
Eigen::VectorXd SurfaceSmoothing::curvatures(ptrdiff_t(0));
int ID = 0, iteration = 0;
int no_boundary_vertices, no_boundary_adjacent_vertices;
double SurfaceSmoothing::vertex_weight = 10.0;
double SurfaceSmoothing::edge_weight = 1.0;
double SurfaceSmoothing::factor = 1.9;
vector<vector<int>> neighbors;


//Map to the mesh IDs
std::unordered_map<int, Eigen::MatrixXd> SurfaceSmoothing::precomputed_L; 
std::unordered_map<int, Eigen::SparseMatrix<double>> SurfaceSmoothing::precompute_matrix_for_LM_and_edges;
std::unordered_map<int, Eigen::SparseMatrix<double>> SurfaceSmoothing::AT_for_LM_and_edges;
std::unordered_map<int, Eigen::SparseMatrix<double>> SurfaceSmoothing::precompute_matrix_for_positions;
std::unordered_map<int, Eigen::SparseMatrix<double>> SurfaceSmoothing::AT_for_positions;
std::unordered_map<int, Eigen::VectorXi> SurfaceSmoothing::precomputed_laplacian_weights;

void SurfaceSmoothing::smooth(Mesh& base_mesh, bool& BOUNDARY_IS_DIRTY){
	if (base_mesh.patches.empty()) {
		base_mesh.patches = Patch::init_patches(base_mesh);
	}

	if(prev_vertex_count != base_mesh.V.rows()) { //The mesh topology has changed, so we need to reset the precomputed matrices. If only boundary constraints got added, we can reuse part of the matrices so don't clear it all out but instead overwrite certain parts
        ID++;
		clear_precomputed_matrices();
		prev_vertex_count = base_mesh.V.rows();
		iteration = 0;
		for (int i = 0; i < base_mesh.patches.size(); i++) {
			solver1_array.push_back(new Eigen::SparseLU<Eigen::SparseMatrix<double>>);
			solver2_array.push_back(new Eigen::SparseLU<Eigen::SparseMatrix<double>>);
		}
	}
	if(BOUNDARY_IS_DIRTY) {
		iteration = 0;
	}



	std::cout << "Check that we odn't have duplicate mesh IDs" << std::endl;
	for (int i = 0; i < base_mesh.patches.size(); i++) {
		Patch* patch = (base_mesh.patches[i]);
		std::cout << (*patch).mesh.ID << std::endl;

		smooth_main((*patch).mesh, BOUNDARY_IS_DIRTY);
		//if (base_mesh.patches.size() == 1) {
			(*patch).update_parent_vertex_positions(base_mesh.V);
	//	}
		prev_mesh_ID = (*patch).mesh.ID;
	}

	BOUNDARY_IS_DIRTY = false;
	iteration++;
}

void SurfaceSmoothing::clear_precomputed_matrices() {
	precomputed_L.clear(); //This reduces the size of the unordered_map to 0

	for(auto it = precompute_matrix_for_LM_and_edges.begin(); it != precompute_matrix_for_LM_and_edges.end(); ++it) {
		it->second.resize(0, 0); //Resize the referenced Eigen::MatrixXd to size 0
	}
	precompute_matrix_for_LM_and_edges.clear();

	for(auto it = precompute_matrix_for_positions.begin(); it != precompute_matrix_for_positions.end(); ++it) {
		it->second.resize(0, 0); //Resize the referenced Eigen::MatrixXd to size 0
	}
	precompute_matrix_for_positions.clear(); 
}

void SurfaceSmoothing::smooth_main(Mesh &m, bool BOUNDARY_IS_DIRTY) {
	Eigen::MatrixXd L = get_precomputed_L(m);
	if(L.rows()==0 && L.cols()==0) { //Stuff that gets computed once, when we get a new mesh topology
		L = compute_laplacian_matrix(m);
		set_precomputed_L(m, L);
		adjacency_list(m.F, neighbors);
	}
	else if (m.ID != prev_mesh_ID){
		adjacency_list(m.F, neighbors);
	}

	if(BOUNDARY_IS_DIRTY || m.ID != prev_mesh_ID) {
		no_boundary_vertices = (m.vertex_boundary_markers.array() > 0).count();
		no_boundary_adjacent_vertices = 0;
		for(int i = 0; i < m.V.rows(); i++) {
			if(m.vertex_boundary_markers[i] > 0) {
				for(int j = 0; j < neighbors[i].size(); j++) {
					if(m.vertex_boundary_markers[neighbors[i][j]] == 0) {
						no_boundary_adjacent_vertices++;
					}
				}
			}
		}
	}
	igl::per_vertex_normals(m.V, m.F, PER_VERTEX_NORMALS_WEIGHTING_TYPE_UNIFORM, vertex_normals);
	Eigen::VectorXd target_LMs = compute_target_LMs(m, L, BOUNDARY_IS_DIRTY);
	Eigen::VectorXd target_edge_lengths = compute_target_edge_lengths(m, L);
	compute_target_vertices(m, L, target_LMs, target_edge_lengths, BOUNDARY_IS_DIRTY);
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
	Eigen::VectorXd initial_curvatures(m.V.rows());
	initial_curvatures.setZero();
	int original_stroke_no_vertices = (m.part_of_original_stroke.array() > 0).count();

	for(int i = 0; i < m.V.rows(); i++) {
		if((m.part_of_original_stroke[i]) && (m.vertex_boundary_markers[i] > 0)) { //Boundary vertex on original stroke
			Eigen::RowVector3d vec(0, 0, 0);
			int count = 0;
			for (int j = 0; j < neighbors[i].size(); j++) {
				if (m.vertex_boundary_markers[neighbors[i][j]] == m.vertex_boundary_markers[i]) {//the neighbor is a boundary vertex on the same stroke
					vec += m.V.row(neighbors[i][j]);
					count++;
				}
			}
			initial_curvatures[i] = (m.V.row(i) - (1.0 / count)*vec).norm(); //TEST: SMOOTH AFTER PULL. maybe revert this change?
			//initial_curvatures[i] = (m.V.row(i) - 0.5 * (m.V.row((i - 1 + original_stroke_no_vertices) % original_stroke_no_vertices) + m.V.row((i + 1 + original_stroke_no_vertices) % original_stroke_no_vertices))).norm(); //Take modulo of the number of original stroke vertices to avoid wrapping the first and last of the original stroke to the interior of the mesh
		} else if(m.vertex_boundary_markers[i] > 0) { //Boundary vertex on an added stroke
			Eigen::RowVector3d vec(0,0,0);
			int count = 0;
			for(int j = 0; j < neighbors[i].size(); j++) {
				if(m.vertex_boundary_markers[neighbors[i][j]] == m.vertex_boundary_markers[i]) {//the neighbor is a boundary vertex on the same stroke
					vec += m.V.row(neighbors[i][j]);
					count++;
				}
			}
			initial_curvatures[i] = (m.V.row(i) - (1.0 / count)*vec).norm();
		}
	}

	return initial_curvatures;
}

Eigen::VectorXd SurfaceSmoothing::compute_target_LMs(Mesh &m, Eigen::MatrixXd &L, bool BOUNDARY_IS_DIRTY) {
	Eigen::SparseMatrix<double> A = get_precompute_matrix_for_LM_and_edges(m);
    Eigen::SparseMatrix<double> AT;
	if(A.rows() == 0 && A.cols() == 0) {
		std::vector<T> tripletList;
		tripletList.reserve(m.V.rows()*(m.V.rows() + 1));
		for(int i = 0; i < m.V.rows(); i++) {
			for (int j = 0; j < m.V.rows(); j++) {
				tripletList.push_back(T(i, j, L(i, j)));
			}
		
			if(m.vertex_boundary_markers[i] > 0) { //Constrain only the boundary in the first iteration
				tripletList.push_back(T(m.V.rows() + i, i, 1));
			}
		}
		A = Eigen::SparseMatrix<double>(m.V.rows() * 2, m.V.rows());
		A.setFromTriplets(tripletList.begin(), tripletList.end());
		A.prune(0.0);
		AT = A.transpose();
		(*solver1_array[m.ID - 1]).compute(AT*A); //Take ID-1 because ID 0 is reserved for the basemesh and patches start from 1
		//solver1.compute(AT*A);
		set_precompute_matrix_for_LM_and_edges(m, A);
        set_AT_for_LM_and_edges(m, AT);
	}
	else if(BOUNDARY_IS_DIRTY && iteration == 0) {
		int nr_V = m.V.rows();
		A.prune([nr_V](int i, int j, float) {return i < nr_V; }); //Keeps the first #V rows (which contains L, which didn't change) and prunes the rest
	
		for(int i = 0; i < m.V.rows(); i++) {
			if(m.vertex_boundary_markers[i] > 0) { //Constrain only the current boundary in the first iteration
				A.coeffRef(m.V.rows() + i, i) = 1; //For target LM'/edge length'
			}
		}

		AT = A.transpose();
		(*solver1_array[m.ID - 1]).compute(AT*A); //Take ID-1 because ID 0 is reserved for the basemesh and patches start from 1
		//solver1.compute(AT*A);
		set_precompute_matrix_for_LM_and_edges(m, A);
		set_AT_for_LM_and_edges(m, AT);
	}
	
	if(iteration == 1) { //Constrain all vertices after the first iteration
		std::vector<T> tripletList;
		tripletList.reserve(m.V.rows()*(m.V.rows() + 1));
		for (int i = 0; i < m.V.rows(); i++) {
			for (int j = 0; j < m.V.rows(); j++) {
				tripletList.push_back(T(i, j, L(i, j)));
			}
			tripletList.push_back(T(m.V.rows() + i, i, 1));
		}

		A.setFromTriplets(tripletList.begin(), tripletList.end());
		A.prune(0.0);
		AT = A.transpose();
		(*solver1_array[m.ID - 1]).compute(AT*A); //Take ID-1 because ID 0 is reserved for the basemesh and patches start from 1
		//solver1.compute(AT*A);
		set_precompute_matrix_for_LM_and_edges(m, A);
        set_AT_for_LM_and_edges(m, AT);
	}
	
    Eigen::VectorXd current_curvatures;
	if(iteration == 0) {
		initial_curvature = compute_initial_curvature(m);
    }else{
        current_curvatures = get_curvatures(m);
    }

	Eigen::VectorXd b = Eigen::VectorXd::Zero(m.V.rows()*2);
    for(int i = 0; i < m.V.rows(); i++) {
		if(iteration == 0) {
			if(m.vertex_boundary_markers[i] > 0) {
				b[m.V.rows() + i] = initial_curvature[i];
			}
		} else {
			b[m.V.rows() + i] = current_curvatures[i];
		}
	}

	AT = get_AT_for_LM_and_edges(m);
	Eigen::VectorXd target_LM = (*solver1_array[m.ID - 1]).solve(AT*b);
	//Eigen::VectorXd target_LM = solver1.solve(AT*b);

	return target_LM;
}

Eigen::VectorXd SurfaceSmoothing::compute_target_edge_lengths(Mesh &m, Eigen::MatrixXd &L) {
	Eigen::SparseMatrix<double> A = get_precompute_matrix_for_LM_and_edges(m);
    Eigen::SparseMatrix<double> AT = get_AT_for_LM_and_edges(m);
	//A will be set up when computing the target LMs

	Eigen::VectorXd b = Eigen::VectorXd::Zero(m.V.rows() * 2);

	double cum_edge_length;
	if(iteration == 0){
		for(int i = 0; i < m.V.rows(); i++) {
			if(m.vertex_boundary_markers[i] > 0) {
				cum_edge_length = 0.0;
				for(int j = 0; j < neighbors[i].size(); j++) {
                    if(m.vertex_boundary_markers[neighbors[i][j]] > 0){
                        cum_edge_length += (m.V.row(neighbors[i][j]) - m.V.row(i)).norm();
                    }
				}
				b[m.V.rows() + i] = (double)(cum_edge_length / neighbors[i].size());
			}
		}
	} else {
		for(int i = 0; i < m.V.rows(); i++) {
			cum_edge_length = 0.0;
			for(int j = 0; j < neighbors[i].size(); j++) {
				cum_edge_length += (m.V.row(neighbors[i][j]) - m.V.row(i)).norm();
			}
			b[m.V.rows() + i] = (double)(cum_edge_length / neighbors[i].size());
		}
	}

	Eigen::VectorXd target_edge_lengths = (*solver1_array[m.ID - 1]).solve(AT*b);
	//Eigen::VectorXd target_edge_lengths = solver1.solve(AT*b);
	return target_edge_lengths;
}

void SurfaceSmoothing::compute_target_vertices(Mesh &m, Eigen::MatrixXd &L, Eigen::VectorXd &target_LMs, Eigen::VectorXd &target_edge_lengths, bool BOUNDARY_IS_DIRTY) {
	Eigen::SparseMatrix<double> A = get_precompute_matrix_for_positions(m);
    Eigen::SparseMatrix<double> AT = get_AT_for_positions(m);
	Eigen::VectorXi laplacian_weights = get_precomputed_laplacian_weights(m);
	if((A.rows() == 0 && A.cols() == 0) || BOUNDARY_IS_DIRTY) {
		laplacian_weights = Eigen::VectorXi::Ones(m.V.rows());
		Eigen::VectorXi on_border = points_on_border(m);
		for (int i = 0; i < on_border.rows(); i++) {
			if (on_border[i]) {
				laplacian_weights[i] *= factor;
			}
		}
		set_precomputed_laplacian_weights(m, laplacian_weights);

		std::vector<T> tripletList;
		tripletList.reserve(m.V.rows()*(m.V.rows() + 13));
		int count = 0, count2 = 0;

		for(int i = 0; i < m.V.rows(); i++) {
			for(int j = 0; j < m.V.rows(); j++) {
				tripletList.push_back(T(i, j, L(i, j)*laplacian_weights[i]));
			}
			if(m.vertex_boundary_markers[i] > 0) { //Constrain only the boundary LMs and edges
				tripletList.push_back(T(m.V.rows() + count, i, vertex_weight));
				count++;

				for(int j = 0; j < neighbors[i].size(); j++) {
					if(m.vertex_boundary_markers[neighbors[i][j]] == 0) {
						tripletList.push_back(T(m.V.rows() + no_boundary_vertices + count2, i, edge_weight));
						tripletList.push_back(T(m.V.rows() + no_boundary_vertices + count2, neighbors[i][j], -edge_weight));
						count2++;
					}
				}
			}
		}

		A = Eigen::SparseMatrix<double>(m.V.rows() + no_boundary_vertices + no_boundary_adjacent_vertices, m.V.rows());
		A.setFromTriplets(tripletList.begin(), tripletList.end());
		A.prune(0.0);
		AT = A.transpose();
		(*solver2_array[m.ID - 1]).compute(AT*A);
		//solver2.compute(AT*A);
		set_precompute_matrix_for_positions(m, A);
        set_AT_for_positions(m, AT);
	}
    
    Eigen::VectorXd bx = Eigen::VectorXd::Zero(m.V.rows() + no_boundary_vertices + no_boundary_adjacent_vertices);
	Eigen::VectorXd by = Eigen::VectorXd::Zero(m.V.rows() + no_boundary_vertices + no_boundary_adjacent_vertices);
	Eigen::VectorXd bz = Eigen::VectorXd::Zero(m.V.rows() + no_boundary_vertices + no_boundary_adjacent_vertices);

	int count = 0, count2 = 0;
	Eigen::Vector3d delta;
	for(int i = 0; i < m.V.rows(); i++) {
        delta = target_LMs[i] * vertex_normals.row(i);

		bx[i] = delta(0) * laplacian_weights[i];
		by[i] = delta(1) * laplacian_weights[i];
		bz[i] = delta(2) * laplacian_weights[i];

		if(m.vertex_boundary_markers[i] > 0) {
			bx[m.V.rows() + count2] = m.V(i, 0) * vertex_weight;
			by[m.V.rows() + count2] = m.V(i, 1) * vertex_weight;
			bz[m.V.rows() + count2] = m.V(i, 2) * vertex_weight;
			count2++;

			for(int j = 0; j < neighbors[i].size(); j++) {
				if(m.vertex_boundary_markers[neighbors[i][j]] == 0) {
					Eigen::Vector3d vec = (m.V.row(i) - m.V.row(neighbors[i][j]));
					vec.normalize();
					Eigen::Vector3d edge_vec = (target_edge_lengths[i] + target_edge_lengths[neighbors[i][j]]) / 2.0 * vec;

					bx[m.V.rows() + no_boundary_vertices + count] = edge_vec(0) * edge_weight;
					by[m.V.rows() + no_boundary_vertices + count] = edge_vec(1) * edge_weight;
					bz[m.V.rows() + no_boundary_vertices + count] = edge_vec(2) * edge_weight;
					count++;
				}
			}
		}
	}

	/*Eigen::VectorXd Vnewx = solver2.solve(AT*bx);
	Eigen::VectorXd Vnewy = solver2.solve(AT*by);
	Eigen::VectorXd Vnewz = solver2.solve(AT*bz);*/
	Eigen::VectorXd Vnewx = (*solver2_array[m.ID - 1]).solve(AT*bx);
	Eigen::VectorXd Vnewy = (*solver2_array[m.ID - 1]).solve(AT*by);
	Eigen::VectorXd Vnewz = (*solver2_array[m.ID - 1]).solve(AT*bz);

	for(int i = 0; i < m.V.rows(); i++) {
		if(m.vertex_boundary_markers[i]==0){ //Only update non-fixed points (curve points are fixed)
			m.V.row(i) << Vnewx[i], Vnewy[i], Vnewz[i];
		}
		else {
			std::cout << "Vert " << i << " is fixed" << std::endl;
		}
	}
}

Eigen::MatrixX3d SurfaceSmoothing::compute_vertex_laplacians(Mesh &m) {
	Eigen::MatrixX3d laplacians(m.V.rows(),3);
    Eigen::Vector3d vec;
    int nr_neighbors;
    for(int i = 0; i < m.V.rows(); i++) {
		vec = Eigen::Vector3d::Zero();
        nr_neighbors = neighbors[i].size();
		for(int j = 0; j < nr_neighbors; j++) {
			vec += m.V.row(neighbors[i][j]);
		}
		laplacians.row(i) = nr_neighbors * m.V.row(i) - vec.transpose();
	}
	return laplacians;
}

Eigen::VectorXd SurfaceSmoothing::get_curvatures(Mesh &m) {
	Eigen::VectorXd curvatures(m.V.rows());
	Eigen::MatrixX3d laplacians = compute_vertex_laplacians(m);
    int sign;
	for(int i = 0; i < m.V.rows(); i++) {
        sign = (laplacians.row(i).dot(vertex_normals.row(i).transpose()) > 0) ? 1 : -1;
        curvatures[i] = laplacians.row(i).norm()*sign;
	}
	return curvatures;
}


//TODO: instead iterate over sharp_edges, and lookup in EV the vertices that belong to these sharp edges. Set point_on_border to true for these
Eigen::VectorXi SurfaceSmoothing::points_on_border(Mesh& m) {
    Eigen::MatrixXi EV, FE, EF;
    igl::edge_topology(m.V, m.F, EV, FE, EF);

	int equal_pos;
	Eigen::VectorXi col1Equals, col2Equals;
	Eigen::VectorXi point_on_border(neighbors.size());
	point_on_border.setZero();

	for(int i = 0; i < neighbors.size(); i++) {
		if(point_on_border[i] == 1) {
			continue;
		}

		for(int j = 0; j < neighbors[i].size(); j++) {
			col1Equals = EV.col(0).cwiseEqual(std::min(i, neighbors[i][j])).cast<int>();
			col2Equals = EV.col(1).cwiseEqual(std::max(i, neighbors[i][j])).cast<int>();
			int maxval = (col1Equals + col2Equals).maxCoeff(&equal_pos); //Find the row that contains both vertices of this edge

			if(m.sharp_edge[equal_pos]) {
				point_on_border[i] = 1;
				point_on_border[neighbors[i][j]] = 1;
				break;
			}
		}
	}

	return point_on_border;
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

Eigen::SparseMatrix<double> SurfaceSmoothing::get_AT_for_LM_and_edges(Mesh &m) {
	return AT_for_LM_and_edges[m.ID];
}

void SurfaceSmoothing::set_AT_for_LM_and_edges(Mesh &m, Eigen::SparseMatrix<double> AT) {
	AT_for_LM_and_edges[m.ID] = AT;
}

Eigen::SparseMatrix<double> SurfaceSmoothing::get_precompute_matrix_for_positions(Mesh &m) {
	return precompute_matrix_for_positions[m.ID];
}

void SurfaceSmoothing::set_precompute_matrix_for_positions(Mesh &m, Eigen::SparseMatrix<double> pos) {
	precompute_matrix_for_positions[m.ID] = pos;
}

Eigen::SparseMatrix<double> SurfaceSmoothing::get_AT_for_positions(Mesh &m) {
	return AT_for_positions[m.ID];
}

void SurfaceSmoothing::set_AT_for_positions(Mesh &m, Eigen::SparseMatrix<double> AT) {
	AT_for_positions[m.ID] = AT;
}

Eigen::VectorXi SurfaceSmoothing::get_precomputed_laplacian_weights(Mesh &m) {
	return precomputed_laplacian_weights[m.ID];
}

void SurfaceSmoothing::set_precomputed_laplacian_weights(Mesh &m, Eigen::VectorXi &laplacian_weights_val) {
	precomputed_laplacian_weights[m.ID] = laplacian_weights_val;
}