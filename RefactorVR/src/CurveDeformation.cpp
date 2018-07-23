#include "CurveDeformation.h"
#include <igl/adjacency_list.h>
#include <igl/edge_topology.h>
#include <Eigen/Sparse>

using namespace std;
using namespace igl;

typedef Eigen::Triplet<double> T;

int moving_vertex_ID, CONSTRAINT_WEIGHT = 10000;// , stroke_ID, no_vertices;// , no_ROI_vert = -1;
double current_ROI_size, curve_diag_length;
//bool CurveDeformation::smooth_deform_mode, stroke_is_loop, prev_loop_type;

vector<Eigen::Matrix3d> Rot;
//vector<int> vert_bindings;

Eigen::RowVector3d start_pos;
//Eigen::VectorXi fixed_indices, fixed_indices_local, is_fixed; //Indicates for every vertex on the pulled curve whether or not it is fixed (aka it is 0 when it is in the ROI and 1 when it is not)
Eigen::VectorXd B, PosRot;
Eigen::MatrixXd original_L0, original_L1;

Eigen::SparseMatrix<double> A, A_L1_T;
Eigen::SparseLU<Eigen::SparseMatrix<double>> solverL1; //Solver for final vertex positions (with L1)
Eigen::SparseLU<Eigen::SparseMatrix<double>> solverPosRot;

vector<vector<int>> CurveDeformation::neighbors;
Eigen::MatrixXi CurveDeformation::EV, CurveDeformation::FE, CurveDeformation::EF;

int prev_range_size = -1;
double DRAG_SCALE = 1.5;
Eigen::VectorXi visited;
Eigen::VectorXd distance_to_vert;
vector<CurveDeformation::PulledCurve> curves;
vector<vector<int>> vertex_triplet_to_rot_idx;
Eigen::VectorXi is_fixed;
Eigen::MatrixXd original_positions;

std::vector<int> fixed_vertices_local;
std::unordered_map<int, int> vertex_global_to_local;
std::unordered_map<int, int> edge_global_to_local;

void CurveDeformation::startPullCurve(Stroke& _stroke, int _moving_vertex_ID, Eigen::MatrixXd& V, Eigen::MatrixXi& F) {
	moving_vertex_ID = _moving_vertex_ID;
	//start_pos = (_stroke.get3DPoints()).row(_handle_ID);
	//moving_vertex_ID = _stroke.get_vertex_idx_for_point(_handle_ID); //The global vertex index (in V) for the moving vertex
	start_pos = V.row(moving_vertex_ID);
	//no_vertices = _stroke.get3DPoints().rows() - 1; //We should ignore the last one, which is a copy of the first one
//	vert_bindings = _stroke.get_closest_vert_bindings();
	//stroke_is_loop = _stroke.is_loop;
	//stroke_ID = _stroke.get_ID();
	curve_diag_length = compute_curve_diag_length(_stroke);
	//handle_ID = _handle_ID;
	//Rot.resize(no_vertices);
	//no_ROI_vert = -1;
	prev_range_size = -1;
	current_ROI_size = 0.0;
	original_positions = V;
	adjacency_list(F, neighbors);
	igl::edge_topology(V, F, EV, FE, EF);
	visited.resize(V.rows());
	distance_to_vert.resize(V.rows());
}


double CurveDeformation::compute_curve_diag_length(Stroke& _stroke) {
	Eigen::Vector3d maxBB = _stroke.get3DPoints().colwise().maxCoeff();
	Eigen::Vector3d minBB = _stroke.get3DPoints().colwise().minCoeff();
	return (maxBB - minBB).norm();
}

void CurveDeformation::pullCurveTest(const Eigen::RowVector3d& pos, Eigen::MatrixXd& V, Eigen::VectorXi& edge_boundary_markers) {
	double drag_size = (pos - start_pos).norm();
	drag_size *= DRAG_SCALE;
	//drag_size /= curve_diag_length; //TODO: Not sure if this is neccessary? pUllmesh2 doesn't do it
	bool ROI_is_updated = false;
	if (!current_ROI_size || (fabs(drag_size - current_ROI_size) > 0.01)) { //Take the current drag_size and current_ROI_size relative to the size of the stroke we're pulling on. //TODO: test if we want to take it relative to the size of the whole mesh instead (we have access to V after all)?
		std::cout << "updating ROI?" << std::endl;
		ROI_is_updated = update_ROI_test(drag_size, V, edge_boundary_markers);
	}

	if (prev_range_size <= 1) {
		V.row(moving_vertex_ID) = pos;
	}
	else {
		if (ROI_is_updated) {
			for (int i = 0; i < curves.size(); i++) {
				setup_for_update_curve_test(curves[i], V);
			}
		}
		V.row(moving_vertex_ID) = pos;
		for (int i = 0; i < curves.size(); i++) {
			update_curve_test(curves[i], V);
		}
	}
}



bool CurveDeformation::update_ROI_test(double drag_size, Eigen::MatrixXd& V, Eigen::VectorXi& edge_boundary_markers) {
	current_ROI_size = drag_size;


	//TODO: PullMesh2 returns all vertices to the positions that they had at the time of being added to the last ROI (so e.g. the handle will have a new "original" position everytime new vertices get added to the ROI
	vector<int> vertices_in_range = collect_vertices_within_drag_length(drag_size, V, edge_boundary_markers);
	if (vertices_in_range.size() == prev_range_size) {
		std::cout << "No change in range size" << std::endl;
		return false;
	}
	prev_range_size = vertices_in_range.size();
	std::cout << "vertices in range: " << prev_range_size << std::endl;
	sort_by_distance(vertices_in_range);

	//Collect vertices in each each curve and construct cuve hierarchy
	Eigen::VectorXi edge_consumed(EV.rows()), fixed(V.rows());
	edge_consumed.setZero();
	fixed.setZero();
	curves.clear();
	for (int i = 0; i < vertices_in_range.size(); i++) {
		while (true) {
			PulledCurve curve;
			bool success = create_pulled_curve_by_propagation(vertices_in_range[i], edge_consumed, curve, edge_boundary_markers);
			if (!success) {
				break;
			}
			curves.push_back(curve);

			for (int j = 0; j < curve.vertices.size(); j++) {
				if (fixed[curve.vertices[j]] && std::find(curve.fixed_vertices.begin(), curve.fixed_vertices.end(), curve.vertices[j]) == curve.fixed_vertices.end()) {
					curve.fixed_vertices.push_back(curve.vertices[j]);
				}
				fixed[curve.vertices[j]] = 1;
			}
		}
	}
	return true;
}

std::vector<int> CurveDeformation::collect_vertices_within_drag_length(double drag_size, Eigen::MatrixXd& V, Eigen::VectorXi& edge_boundary_markers) {
	visited.setZero();
	distance_to_vert.setZero();
	visited[moving_vertex_ID] = 1;
	distance_to_vert[moving_vertex_ID] = drag_size;
	std::cout << "moving vertex id" << moving_vertex_ID << std::endl;
	vector<int> vertices_in_range;
	vertices_in_range.push_back(moving_vertex_ID);
	propagate(-1, moving_vertex_ID, drag_size, vertices_in_range, V, edge_boundary_markers);
	std::cout << "Vertices in range: ";
	for (int i = 0; i < vertices_in_range.size(); i++) {
		std::cout << vertices_in_range[i] << "   ";
	}
	return vertices_in_range;
}

void CurveDeformation::propagate(int prev_edge, int vert, double remaining_distance, std::vector<int>& vertices_in_range, Eigen::MatrixXd& V, Eigen::VectorXi& edge_boundary_markers) {
	int next_vert;
	for (int i = 0; i < neighbors[vert].size(); i++) {
		int edge = find_edge(vert, neighbors[vert][i]);
		if (edge == prev_edge) {
			continue;
		}
		double edge_length = (original_positions.row(vert) - original_positions.row(neighbors[vert][i])).norm();
		//double edge_length = (V.row(vert) - V.row(neighbors[vert][i])).norm();
		if (edge_boundary_markers[edge]) {
			if (remaining_distance - edge_length > 0) {
				next_vert = neighbors[vert][i];
				if (visited[next_vert]) {
					continue;
				}
				visited[next_vert] = 1;
				distance_to_vert[next_vert] = remaining_distance - edge_length;
				vertices_in_range.push_back(next_vert);
				propagate(edge, next_vert, remaining_distance - edge_length, vertices_in_range, V, edge_boundary_markers);
			}
		}
	}
}

int CurveDeformation::find_edge(int start, int end) {
	Eigen::VectorXi col1Equals, col2Equals;
	int equal_pos;
	col1Equals = EV.col(0).cwiseEqual(min(start, end)).cast<int>();
	col2Equals = EV.col(1).cwiseEqual(max(start, end)).cast<int>();
	(col1Equals + col2Equals).maxCoeff(&equal_pos); //Find the row that contains both vertices of this edge

	return equal_pos;
}

void CurveDeformation::sort_by_distance(std::vector<int>& vertices_in_range) {
	sort(vertices_in_range.begin(), vertices_in_range.end(), distance_sorter);
}

bool CurveDeformation::distance_sorter(int a, int b) {
	return distance_to_vert[a] > distance_to_vert[b];
}

bool CurveDeformation::create_pulled_curve_by_propagation(int vert, Eigen::VectorXi& edge_consumed, PulledCurve& curve, Eigen::VectorXi& edge_boundary_markers) {
	curve.vertices.push_back(vert);
	curve.fixed_vertices.push_back(vert);

	int edge_marker = 0;
	for (int i = 0; i < neighbors[vert].size(); i++) {
		int edge = find_edge(vert, neighbors[vert][i]);
		if (edge_boundary_markers[edge] && !edge_consumed[edge] && visited[neighbors[vert][i]]) {
			edge_marker = edge_boundary_markers[edge];
		}
	}

	if (!edge_marker) {
		return false;
	}

	propagate(curve, -1, vert, edge_marker, edge_consumed, edge_boundary_markers);

	for (int i = 0; i < curve.edges.size(); i++) {
		//TODO: check if I need this
		int prev_edge = get_adjacent_seam_edge(EV(curve.edges[i], 0), curve.edges[i], curve.edges, edge_boundary_markers);
		int next_edge = get_adjacent_seam_edge(EV(curve.edges[i], 1), curve.edges[i], curve.edges, edge_boundary_markers);
		if (prev_edge != -1 && next_edge != -1) {
			curve.edge_triplets.conservativeResize(curve.edge_triplets.rows() + 1, Eigen::NoChange);
			curve.edge_triplets.bottomRows(1) << prev_edge, curve.edges[i], next_edge;
		}
	}

	for (int i = 0; i < curve.vertices.size(); i++) {
		vector<int> adjacent_vertices = get_adjacent_seam_vertices(curve.vertices[i], curve);
		if (adjacent_vertices.size() == 2) {
			curve.vertex_triplets.conservativeResize(curve.vertex_triplets.rows() + 1, Eigen::NoChange);
			curve.vertex_triplets.bottomRows(1) << adjacent_vertices[0], curve.vertices[i], adjacent_vertices[1];
		}
	}
	return true;
}

vector<int> CurveDeformation::get_adjacent_seam_vertices(int vert, PulledCurve& curve) {
	vector<int> adjacent_vertices;
	for (int i = 0; i < neighbors[vert].size(); i++) {
		int edge = find_edge(vert, neighbors[vert][i]);
		if (std::find(curve.edges.begin(), curve.edges.end(), edge) != curve.edges.end()) {
			adjacent_vertices.push_back(EV(edge, 0) == vert ? EV(edge, 1) : EV(edge, 0));
			if (adjacent_vertices.size() == 2) {
				return adjacent_vertices;
			}
		}
	}
	return adjacent_vertices;
}

int CurveDeformation::get_adjacent_seam_edge(int vert, int edge, std::vector<int> edges, Eigen::VectorXi & edge_boundary_markers) {
	for (int i = 0; i < neighbors[vert].size(); i++) {
		int edge2 = find_edge(vert, neighbors[vert][i]);
		if (edge2 != edge && edge_boundary_markers[edge] == edge_boundary_markers[edge2]) {
			if (std::find(edges.begin(), edges.end(), edge2) != edges.end()) {
				return edge2;
			}
		}
	}
	return -1;
}

void CurveDeformation::propagate(PulledCurve& curve, int prev_edge, int vert, int edge_marker, Eigen::VectorXi& edge_consumed, Eigen::VectorXi& edge_boundary_markers) {
	bool is_terminal = true;
	for (int i = 0; i < neighbors[vert].size(); i++) {
		int edge = find_edge(vert, neighbors[vert][i]);
		if (edge != prev_edge && edge_boundary_markers[edge] == edge_marker && !edge_consumed[edge]) {
			if (visited[neighbors[vert][i]]) {
				curve.edges.push_back(edge);
				edge_consumed[edge] = 1;
				if (std::find(curve.vertices.begin(), curve.vertices.end(), neighbors[vert][i]) != curve.vertices.end()) { //Curve vertices already contain the next vertex, so we've gotten a loop
					return;
				}

				curve.vertices.push_back(neighbors[vert][i]);
				propagate(curve, edge, neighbors[vert][i], edge_marker, edge_consumed, edge_boundary_markers);

				is_terminal = false;
			}
		}
	}

	if (is_terminal) {
		curve.fixed_vertices.push_back(vert);
		if (prev_edge == -1) {
			return;
		}
		int next_edge = find_next_edge(vert, prev_edge, edge_marker, edge_boundary_markers);
		if (next_edge == -1) {
			return;
		}

		if (!edge_consumed[next_edge]) {
			curve.edges.push_back(next_edge);
			edge_consumed[next_edge] = 1;
		}
		if (std::find(curve.fixed_edges.begin(), curve.fixed_edges.end(), next_edge) == curve.fixed_edges.end()) {
			curve.fixed_edges.push_back(next_edge);
		}

		int next_vert = EV(next_edge, 0) == vert ? EV(next_edge, 1) : EV(next_edge, 0); //Get the vertex on the other side of the edge
		if (std::find(curve.vertices.begin(), curve.vertices.end(), next_vert) == curve.vertices.end()) {
			curve.vertices.push_back(next_vert);
		}
		if (std::find(curve.fixed_vertices.begin(), curve.fixed_vertices.end(), next_vert) == curve.fixed_vertices.end()) {
			curve.fixed_vertices.push_back(next_vert);
		}
	}
}

int CurveDeformation::find_next_edge(int vert, int prev_edge, int edge_marker, Eigen::VectorXi& edge_boundary_markers) {
	int edge;
	for (int i = 0; i < neighbors[vert].size(); i++) {
		edge = find_edge(vert, neighbors[vert][i]);
		if (edge != prev_edge && edge_marker == edge_boundary_markers[edge]) {
			return edge;
		}
	}
	return -1;
}

void CurveDeformation::setup_for_update_curve_test(PulledCurve& curve, Eigen::MatrixXd& V) {
	A.resize(curve.edges.size() * 3 + curve.edge_triplets.rows() * 9 + curve.fixed_vertices.size() * 3 + curve.fixed_edges.size() * 3, curve.vertices.size() * 3 + curve.edges.size() * 3);
	B.resize(curve.edges.size() * 3 + curve.edge_triplets.rows() * 9 + curve.fixed_vertices.size() * 3 + curve.fixed_edges.size() * 3);
	Rot.resize(curve.edges.size());
	original_L0.resize(curve.edges.size(), 3);
	for (int i = 0; i < curve.edges.size(); i++) {
		original_L0.row(i) = V.row(EV(curve.edges[i],0)) - V.row(EV(curve.edges[i],1));
	}

	for (int i = 0; i < curve.vertices.size(); i++) {
		vertex_global_to_local.insert({ curve.vertices[i], i });
	}

	int idx;
	fixed_vertices_local.resize(curve.fixed_vertices.size());
	for (int i = 0; i < curve.fixed_vertices.size(); i++) { //Create a copy with local indexing for fixed_vertices, and also keep the global indexing
		auto loc = std::find(curve.vertices.begin(), curve.vertices.end(), curve.fixed_vertices[i]);
		idx = distance(curve.vertices.begin(), loc);
		fixed_vertices_local[i] = idx;
	}

	std::cout << curve.edges.size() <<"   Curve fixed edges before: " << std::endl;
	for (int i = 0; i < curve.fixed_edges.size(); i++) { //Create local indexing for fixed_edges
		std::cout << curve.fixed_edges[i] << "  ";
		auto loc = std::find(curve.edges.begin(), curve.edges.end(), curve.fixed_edges[i]);
		idx = distance(curve.edges.begin(), loc);
		curve.fixed_edges[i] = idx;
	}
	std::cout << "After: " << std::endl;
	for (int i = 0; i < curve.fixed_edges.size(); i++) {
		std::cout << curve.fixed_edges[i] << " ";
	}
	cout << std::endl;

	for (int i = 0; i < curve.edges.size(); i++) {
		edge_global_to_local.insert({ curve.edges[i], i });
	}

	is_fixed.resize(curve.vertices.size());
	is_fixed.setZero();
	for (int i = 0; i < fixed_vertices_local.size(); i++) {
		is_fixed[fixed_vertices_local[i]] = 1;
	}

	setup_for_L1_position_step_test(curve, V);
}

void CurveDeformation::setup_for_L1_position_step_test(PulledCurve& curve, Eigen::MatrixXd& V) {
	std::vector<T> tripletList;
	tripletList.reserve(curve.vertex_triplets.rows() * 3 + curve.fixed_vertices.size());
	original_L1.resize(curve.vertex_triplets.rows(), 3);
	vertex_triplet_to_rot_idx.resize(curve.vertex_triplets.rows(), vector<int>(2));

	for (int i = 0; i < curve.vertex_triplets.rows(); i++) {
		int v0 = curve.vertex_triplets(i, 0);
		int v1 = curve.vertex_triplets(i, 1);
		int v2 = curve.vertex_triplets(i, 2);

		original_L1.row(i) = V.row(v1) - (0.5*V.row(v0) + 0.5*V.row(v2));
		vertex_triplet_to_rot_idx[i][0] = edge_global_to_local.find(find_edge(v0, v1))->second;
		vertex_triplet_to_rot_idx[i][1] = edge_global_to_local.find(find_edge(v1, v2))->second;

		tripletList.push_back(T(i, vertex_global_to_local.find(v0)->second, -0.5));
		tripletList.push_back(T(i, vertex_global_to_local.find(v1)->second, 1.0));
		tripletList.push_back(T(i, vertex_global_to_local.find(v2)->second, -0.5));
	}

	for (int i = 0; i < fixed_vertices_local.size(); i++) {
		tripletList.push_back(T(curve.vertex_triplets.rows() + i, fixed_vertices_local[i], CONSTRAINT_WEIGHT));
	}
	Eigen::SparseMatrix<double> A_L1(curve.vertex_triplets.rows() + fixed_vertices_local.size(), curve.vertices.size());
	A_L1.setFromTriplets(tripletList.begin(), tripletList.end());

	A_L1_T = A_L1.transpose();
	solverL1.compute(A_L1_T*A_L1);
}

void CurveDeformation::update_curve_test(PulledCurve& curve, Eigen::MatrixXd& V) {
	for (int i = 0; i < curve.edges.size(); i++) {
		Rot[i] = Eigen::Matrix3d::Identity();
	}

	for (int i = 0; i < 2; i++) {
		solve_for_pos_and_rot_test(V, curve);
		update_rot_test(curve);
	}
	final_L1_pos_test(V, curve);
}

void CurveDeformation::solve_for_pos_and_rot_test(Eigen::MatrixXd& V, PulledCurve& curve) {
	A.setZero();
	B.setZero();
	std::vector<T> tripletList;
	tripletList.reserve(curve.edges.size()*15 + curve.edge_triplets.rows()*3*27 + fixed_vertices_local.size()*3 + curve.fixed_edges.size() * 3);

	int v0, v1;
	for (int i = 0; i < curve.edges.size(); i++) { //Laplacian of the position: v1' - v0' = dR * R * L(v1 - v0)
		v0 = vertex_global_to_local.find(EV(curve.edges[i], 0))->second;
		v1 = vertex_global_to_local.find(EV(curve.edges[i], 1))->second;
		
		tripletList.push_back(T(i * 3 + 0, v0 * 3, -1)); //L0
		tripletList.push_back(T(i * 3 + 0, v1 * 3, 1)); //L0
		tripletList.push_back(T(i * 3 + 0, curve.vertices.size() * 3 + i * 3, 0));
		tripletList.push_back(T(i * 3 + 0, curve.vertices.size() * 3 + i * 3 + 1, -(Rot[i].row(2).dot(original_L0.row(i)))));
		tripletList.push_back(T(i * 3 + 0, curve.vertices.size() * 3 + i * 3 + 2, Rot[i].row(1).dot(original_L0.row(i))));
		B[i * 3] = Rot[i].row(0).dot(original_L0.row(i)); //Constant

		tripletList.push_back(T(i * 3 + 1, v0 * 3 + 1, -1)); //L0
		tripletList.push_back(T(i * 3 + 1, v1 * 3 + 1, 1)); //L0
		tripletList.push_back(T(i * 3 + 1, curve.vertices.size() * 3 + i * 3, Rot[i].row(2).dot(original_L0.row(i))));
		tripletList.push_back(T(i * 3 + 1, curve.vertices.size() * 3 + i * 3 + 1, 0));
		tripletList.push_back(T(i * 3 + 1, curve.vertices.size() * 3 + i * 3 + 2, -Rot[i].row(0).dot(original_L0.row(i))));
		B[i * 3 + 1] = Rot[i].row(1).dot(original_L0.row(i)); //Constant

		tripletList.push_back(T(i * 3 + 2, v0 * 3 + 2, -1)); //L0
		tripletList.push_back(T(i * 3 + 2, v1 * 3 + 2, 1)); //L0
		tripletList.push_back(T(i * 3 + 2, curve.vertices.size() * 3 + i * 3, -Rot[i].row(1).dot(original_L0.row(i)))); //Note: this is different from teddy, but I think teddy is wrong here
		tripletList.push_back(T(i * 3 + 2, curve.vertices.size() * 3 + i * 3 + 1, Rot[i].row(0).dot(original_L0.row(i))));
		tripletList.push_back(T(i * 3 + 2, curve.vertices.size() * 3 + i * 3 + 2, 0));
		B[i * 3 + 2] = Rot[i].row(2).dot(original_L0.row(i)); //Constant
	}

	int e0, e1, e2;
	for (int i = 0; i < curve.edge_triplets.rows(); i++) { //Laplacian of the rotation matrix: L dR*R = 0 (r_i*R_i - r_j*R_j = 0 for i,j in Edges)
		e0 = edge_global_to_local.find(curve.edge_triplets(i, 0))->second;
		e1 = edge_global_to_local.find(curve.edge_triplets(i, 1))->second;
		e2 = edge_global_to_local.find(curve.edge_triplets(i, 2))->second;

		for (int j = 0; j < 3; j++) {
			tripletList.push_back(T(curve.edges.size() * 3 + i * 9 + j, curve.vertices.size() * 3 + e0 * 3 + 0, 0));
			tripletList.push_back(T(curve.edges.size() * 3 + i * 9 + j, curve.vertices.size() * 3 + e0 * 3 + 1, -0.5 * Rot[e0](2, j)));
			tripletList.push_back(T(curve.edges.size() * 3 + i * 9 + j, curve.vertices.size() * 3 + e0 * 3 + 2, -0.5 * -Rot[e0](1, j)));
			tripletList.push_back(T(curve.edges.size() * 3 + i * 9 + j, curve.vertices.size() * 3 + e1 * 3 + 0, 0));
			tripletList.push_back(T(curve.edges.size() * 3 + i * 9 + j, curve.vertices.size() * 3 + e1 * 3 + 1, Rot[e1](2, j)));
			tripletList.push_back(T(curve.edges.size() * 3 + i * 9 + j, curve.vertices.size() * 3 + e1 * 3 + 2, -Rot[e1](1, j)));
			tripletList.push_back(T(curve.edges.size() * 3 + i * 9 + j, curve.vertices.size() * 3 + e2 * 3 + 0, 0));
			tripletList.push_back(T(curve.edges.size() * 3 + i * 9 + j, curve.vertices.size() * 3 + e2 * 3 + 1, -0.5 * Rot[e2](2, j)));
			tripletList.push_back(T(curve.edges.size() * 3 + i * 9 + j, curve.vertices.size() * 3 + e2 * 3 + 2, -0.5 * -Rot[e2](1, j)));
			B[curve.edges.size() * 3 + i * 9 + j] = 0.5 * Rot[e0](0, j) - 1 * Rot[e1](0, j) + 0.5 * Rot[e2](0, j);

			tripletList.push_back(T(curve.edges.size() * 3 + i * 9 + 3 + j, curve.vertices.size() * 3 + e0 * 3 + 0, -0.5 * -Rot[e0](2, j)));
			tripletList.push_back(T(curve.edges.size() * 3 + i * 9 + 3 + j, curve.vertices.size() * 3 + e0 * 3 + 1, 0));
			tripletList.push_back(T(curve.edges.size() * 3 + i * 9 + 3 + j, curve.vertices.size() * 3 + e0 * 3 + 2, -0.5 * Rot[e0](0, j)));
			tripletList.push_back(T(curve.edges.size() * 3 + i * 9 + 3 + j, curve.vertices.size() * 3 + e1 * 3 + 0, -Rot[e1](2, j)));
			tripletList.push_back(T(curve.edges.size() * 3 + i * 9 + 3 + j, curve.vertices.size() * 3 + e1 * 3 + 1, 0));
			tripletList.push_back(T(curve.edges.size() * 3 + i * 9 + 3 + j, curve.vertices.size() * 3 + e1 * 3 + 2, Rot[e1](0, j)));
			tripletList.push_back(T(curve.edges.size() * 3 + i * 9 + 3 + j, curve.vertices.size() * 3 + e2 * 3 + 0, -0.5 * -Rot[e2](2, j)));
			tripletList.push_back(T(curve.edges.size() * 3 + i * 9 + 3 + j, curve.vertices.size() * 3 + e2 * 3 + 1, 0));
			tripletList.push_back(T(curve.edges.size() * 3 + i * 9 + 3 + j, curve.vertices.size() * 3 + e2 * 3 + 2, -0.5 * Rot[e2](0, j)));
			B[curve.edges.size() * 3 + i * 9 + 3 + j] = 0.5 * Rot[e0](1, j) - 1 * Rot[e1](1, j) + 0.5 * Rot[e2](1, j);

			tripletList.push_back(T(curve.edges.size() * 3 + i * 9 + 6 + j, curve.vertices.size() * 3 + e0 * 3 + 0, -0.5*Rot[e0](1,j)));
			tripletList.push_back(T(curve.edges.size() * 3 + i * 9 + 6 + j, curve.vertices.size() * 3 + e0 * 3 + 1, -0.5 * -Rot[e0](0, j)));
			tripletList.push_back(T(curve.edges.size() * 3 + i * 9 + 6 + j, curve.vertices.size() * 3 + e0 * 3 + 2, 0));
			tripletList.push_back(T(curve.edges.size() * 3 + i * 9 + 6 + j, curve.vertices.size() * 3 + e1 * 3 + 0, Rot[e1](1, j)));
			tripletList.push_back(T(curve.edges.size() * 3 + i * 9 + 6 + j, curve.vertices.size() * 3 + e1 * 3 + 1, -Rot[e1](0, j)));
			tripletList.push_back(T(curve.edges.size() * 3 + i * 9 + 6 + j, curve.vertices.size() * 3 + e1 * 3 + 2, 0));
			tripletList.push_back(T(curve.edges.size() * 3 + i * 9 + 6 + j, curve.vertices.size() * 3 + e2 * 3 + 0, -0.5*Rot[e2](1, j)));
			tripletList.push_back(T(curve.edges.size() * 3 + i * 9 + 6 + j, curve.vertices.size() * 3 + e2 * 3 + 1, -0.5 * -Rot[e2](0, j)));
			tripletList.push_back(T(curve.edges.size() * 3 + i * 9 + 6 + j, curve.vertices.size() * 3 + e2 * 3 + 2, 0));
			B[curve.edges.size() * 3 + i * 9 + 6 + j] = 0.5 * Rot[e0](2, j) - 1 * Rot[e1](2, j) + 0.5 * Rot[e2](2, j);
		}
	}

	for (int i = 0; i < fixed_vertices_local.size(); i++) { //Position of fixed vertices: v_i - v_i' = 0
		tripletList.push_back(T(curve.edges.size() * 3 + curve.edge_triplets.rows() * 9 + i * 3, fixed_vertices_local[i] * 3, CONSTRAINT_WEIGHT));
		tripletList.push_back(T(curve.edges.size() * 3 + curve.edge_triplets.rows() * 9 + i * 3 + 1, fixed_vertices_local[i] * 3 + 1, CONSTRAINT_WEIGHT));
		tripletList.push_back(T(curve.edges.size() * 3 + curve.edge_triplets.rows() * 9 + i * 3 + 2, fixed_vertices_local[i] * 3 + 2, CONSTRAINT_WEIGHT));

		B[curve.edges.size() * 3 + curve.edge_triplets.rows() * 9 + i * 3] = CONSTRAINT_WEIGHT * V(curve.fixed_vertices[i], 0);
		B[curve.edges.size() * 3 + curve.edge_triplets.rows() * 9 + i * 3 + 1] = CONSTRAINT_WEIGHT * V(curve.fixed_vertices[i], 1);
		B[curve.edges.size() * 3 + curve.edge_triplets.rows() * 9 + i * 3 + 2] = CONSTRAINT_WEIGHT * V(curve.fixed_vertices[i], 2);
	}

	for (int i = 0; i < curve.fixed_edges.size(); i++) {
		tripletList.push_back(T(curve.edges.size() * 3 + curve.edge_triplets.rows() * 9 + fixed_vertices_local.size() * 3 + i * 3 + 0, curve.vertices.size() * 3 + curve.fixed_edges[i] * 3, CONSTRAINT_WEIGHT));
		tripletList.push_back(T(curve.edges.size() * 3 + curve.edge_triplets.rows() * 9 + fixed_vertices_local.size() * 3 + i * 3 + 1, curve.vertices.size() * 3 + curve.fixed_edges[i] * 3 + 1, CONSTRAINT_WEIGHT));
		tripletList.push_back(T(curve.edges.size() * 3 + curve.edge_triplets.rows() * 9 + fixed_vertices_local.size() * 3 + i * 3 + 2, curve.vertices.size() * 3 + curve.fixed_edges[i] * 3 + 2, CONSTRAINT_WEIGHT));

		B[curve.edges.size() * 3 + curve.edge_triplets.rows() * 9 + fixed_vertices_local.size() * 3 + i * 3 + 0] = 0;
		B[curve.edges.size() * 3 + curve.edge_triplets.rows() * 9 + fixed_vertices_local.size() * 3 + i * 3 + 1] = 0;
		B[curve.edges.size() * 3 + curve.edge_triplets.rows() * 9 + fixed_vertices_local.size() * 3 + i * 3 + 2] = 0;
	}

	A.setFromTriplets(tripletList.begin(), tripletList.end());
	A.prune(0.0);
	Eigen::SparseMatrix<double> AT = A.transpose();
	solverPosRot.compute(AT*A);
	PosRot = solverPosRot.solve(AT*B);
}

/*void CurveDeformation::solve_for_pos_and_rot(Eigen::MatrixXd& V) {
	A.setZero();
	B.setZero();
	std::vector<T> tripletList;
	tripletList.reserve(no_vertices * 9 + no_vertices * 3 * 18 + fixed_indices.size() * 6);

	int prev, cur;
	int start = !stroke_is_loop; //Make sure to skip wrapping around for index 0
	for (int i = 0; i < no_vertices; i++) {
		prev = (((i - 1) + no_vertices) % no_vertices);
		cur = i;
		if (i > 0 || stroke_is_loop) {
			tripletList.push_back(T(i * 3 + 0, prev * 3, -1)); //L0
			tripletList.push_back(T(i * 3 + 0, cur * 3, 1)); //L0
		}
		tripletList.push_back(T(i * 3, no_vertices * 3 + cur * 3, 0));
		tripletList.push_back(T(i * 3, no_vertices * 3 + cur * 3 + 1, -(Rot[i].row(2).dot(original_L0.row(i)))));
		tripletList.push_back(T(i * 3, no_vertices * 3 + cur * 3 + 2, Rot[i].row(1).dot(original_L0.row(i))));
		B[i * 3] = Rot[i].row(0).dot(original_L0.row(i)); //constant

		if (i > 0 || stroke_is_loop) {
			tripletList.push_back(T(i * 3 + 1, prev * 3 + 1, -1)); //L0
			tripletList.push_back(T(i * 3 + 1, cur * 3 + 1, 1)); //L0
		}
		tripletList.push_back(T(i * 3 + 1, no_vertices * 3 + cur * 3, Rot[i].row(2).dot(original_L0.row(i))));
		tripletList.push_back(T(i * 3 + 1, no_vertices * 3 + cur * 3 + 1, 0));
		tripletList.push_back(T(i * 3 + 1, no_vertices * 3 + cur * 3 + 2, -Rot[i].row(0).dot(original_L0.row(i))));
		B[i * 3 + 1] = Rot[i].row(1).dot(original_L0.row(i));

		if (i > 0 || stroke_is_loop) {
			tripletList.push_back(T(i * 3 + 2, prev * 3 + 2, -1)); //L0
			tripletList.push_back(T(i * 3 + 2, cur * 3 + 2, 1)); //L0
		}
		tripletList.push_back(T(i * 3 + 2, no_vertices * 3 + cur * 3, -Rot[i].row(1).dot(original_L0.row(i)))); //Note: this is different from teddy, but I think teddy is wrong here
		tripletList.push_back(T(i * 3 + 2, no_vertices * 3 + cur * 3 + 1, Rot[i].row(0).dot(original_L0.row(i))));
		tripletList.push_back(T(i * 3 + 2, no_vertices * 3 + cur * 3 + 2, 0));
		B[i * 3 + 2] = Rot[i].row(2).dot(original_L0.row(i));
	}

	//Setup for r_i*R_i - r_j*R_j = 0 for i,j in Edges
	int prev_rot;
	for (int i = start; i < no_vertices; i++) { //For non-closed strokes, start at index 1 (we can't fill in a partial entry since they depend on eachother, unlike the previous block)
		prev = (((i - 1) + no_vertices) % no_vertices);
		prev_rot = (((i - 1) + no_vertices) % no_vertices);
		for (int j = 0; j < 3; j++) {
			tripletList.push_back(T(no_vertices * 3 + i * 9 + j, no_vertices * 3 + prev * 3, 0));
			tripletList.push_back(T(no_vertices * 3 + i * 9 + j, no_vertices * 3 + prev * 3 + 1, -1 * Rot[prev_rot](2, j)));
			tripletList.push_back(T(no_vertices * 3 + i * 9 + j, no_vertices * 3 + prev * 3 + 2, -1 * -Rot[prev_rot](1, j)));
			tripletList.push_back(T(no_vertices * 3 + i * 9 + j, no_vertices * 3 + i * 3, 0));
			tripletList.push_back(T(no_vertices * 3 + i * 9 + j, no_vertices * 3 + i * 3 + 1, 1 * Rot[i](2, j)));
			tripletList.push_back(T(no_vertices * 3 + i * 9 + j, no_vertices * 3 + i * 3 + 2, 1 * -Rot[i](1, j)));


			B[no_vertices * 3 + i * 9 + j] = -1 * Rot[i](0, j) + 1 * Rot[prev_rot](0, j);

			tripletList.push_back(T(no_vertices * 3 + i * 9 + 3 + j, no_vertices * 3 + prev * 3, -1 * -Rot[prev_rot](2, j)));
			tripletList.push_back(T(no_vertices * 3 + i * 9 + 3 + j, no_vertices * 3 + prev * 3 + 1, 0));
			tripletList.push_back(T(no_vertices * 3 + i * 9 + 3 + j, no_vertices * 3 + prev * 3 + 2, -1 * Rot[prev_rot](0, j)));
			tripletList.push_back(T(no_vertices * 3 + i * 9 + 3 + j, no_vertices * 3 + i * 3, 1 * -Rot[i](2, j)));
			tripletList.push_back(T(no_vertices * 3 + i * 9 + 3 + j, no_vertices * 3 + i * 3 + 1, 0));
			tripletList.push_back(T(no_vertices * 3 + i * 9 + 3 + j, no_vertices * 3 + i * 3 + 2, 1 * Rot[i](0, j)));

			B[no_vertices * 3 + i * 9 + 3 + j] = -1 * Rot[i](1, j) + 1 * Rot[prev_rot](1, j);

			tripletList.push_back(T(no_vertices * 3 + i * 9 + 6 + j, no_vertices * 3 + prev * 3, -1 * Rot[prev_rot](1, j)));
			tripletList.push_back(T(no_vertices * 3 + i * 9 + 6 + j, no_vertices * 3 + prev * 3 + 1, -1 * -Rot[prev_rot](0, j)));
			tripletList.push_back(T(no_vertices * 3 + i * 9 + 6 + j, no_vertices * 3 + prev * 3 + 2, 0));
			tripletList.push_back(T(no_vertices * 3 + i * 9 + 6 + j, no_vertices * 3 + i * 3, 1 * Rot[i](1, j)));
			tripletList.push_back(T(no_vertices * 3 + i * 9 + 6 + j, no_vertices * 3 + i * 3 + 1, 1 * -Rot[i](0, j)));
			tripletList.push_back(T(no_vertices * 3 + i * 9 + 6 + j, no_vertices * 3 + i * 3 + 2, 0));
			B[no_vertices * 3 + i * 9 + 6 + j] = -1 * Rot[i](2, j) + 1 * Rot[prev_rot](2, j);
		}
	}

	//Setup for v_i - v_i' = 0
	for (int i = 0; i < fixed_indices.size(); i++) {
		tripletList.push_back(T(no_vertices * 3 + no_vertices * 9 + i * 3, fixed_indices_local[i] * 3, CONSTRAINT_WEIGHT));
		tripletList.push_back(T(no_vertices * 3 + no_vertices * 9 + i * 3 + 1, fixed_indices_local[i] * 3 + 1, CONSTRAINT_WEIGHT));
		tripletList.push_back(T(no_vertices * 3 + no_vertices * 9 + i * 3 + 2, fixed_indices_local[i] * 3 + 2, CONSTRAINT_WEIGHT));


		B[no_vertices * 3 + no_vertices * 9 + i * 3] = CONSTRAINT_WEIGHT * V(fixed_indices[i], 0);
		B[no_vertices * 3 + no_vertices * 9 + i * 3 + 1] = CONSTRAINT_WEIGHT * V(fixed_indices[i], 1);
		B[no_vertices * 3 + no_vertices * 9 + i * 3 + 2] = CONSTRAINT_WEIGHT * V(fixed_indices[i], 2);
	}

	//Setup for r_i*R_i - R_i' = 0 (note that we want all off-diagonal elements in r_i to be equal to 0 in order to multiply R_i with the identity matrix)
	for (int i = 0; i < fixed_indices.size() - 1; i++) { //Skip the last fixed index (the handle index)
		tripletList.push_back(T(no_vertices * 3 + no_vertices * 9 + fixed_indices.size() * 3 + i * 3, no_vertices * 3 + fixed_indices_local[i] * 3, CONSTRAINT_WEIGHT));
		tripletList.push_back(T(no_vertices * 3 + no_vertices * 9 + fixed_indices.size() * 3 + i * 3 + 1, no_vertices * 3 + fixed_indices_local[i] * 3 + 1, CONSTRAINT_WEIGHT));
		tripletList.push_back(T(no_vertices * 3 + no_vertices * 9 + fixed_indices.size() * 3 + i * 3 + 2, no_vertices * 3 + fixed_indices_local[i] * 3 + 2, CONSTRAINT_WEIGHT));

		B[no_vertices * 3 + no_vertices * 9 + fixed_indices.size() * 3 + i * 3] = 0;
		B[no_vertices * 3 + no_vertices * 9 + fixed_indices.size() * 3 + i * 3 + 1] = 0;
		B[no_vertices * 3 + no_vertices * 9 + fixed_indices.size() * 3 + i * 3 + 2] = 0;
	}

	A.setFromTriplets(tripletList.begin(), tripletList.end());
	A.prune(0.0);
	Eigen::SparseMatrix<double> AT = A.transpose();
	solverPosRot.compute(AT*A);
	PosRot = solverPosRot.solve(AT*B);
}*/

/*void CurveDeformation::update_rot() {
	Eigen::Matrix3d newRot;
	double rx, ry, rz;
	for (int i = 0; i < Rot.size(); i++) {
		rx = PosRot[no_vertices * 3 + i * 3];
		ry = PosRot[no_vertices * 3 + i * 3 + 1];
		rz = PosRot[no_vertices * 3 + i * 3 + 2];

		newRot.row(0) << 1, -rz, ry;
		newRot.row(1) << rz, 1, -rx;
		newRot.row(2) << -ry, rx, 1;
		Rot[i] = newRot*Rot[i]; //Safe for matrix-matrix multiplication with Eigen
		Rot[i] = compute_orthonormal(Rot[i]);
	}
}*/

void CurveDeformation::update_rot_test(PulledCurve& curve) {
	Eigen::Matrix3d newRot;
	double rx, ry, rz;
	for (int i = 0; i < Rot.size(); i++) {
		rx = PosRot[curve.vertices.size() * 3 + i * 3];
		ry = PosRot[curve.vertices.size() * 3 + i * 3 + 1];
		rz = PosRot[curve.vertices.size() * 3 + i * 3 + 2];

		newRot.row(0) << 1, -rz, ry;
		newRot.row(1) << rz, 1, -rx;
		newRot.row(2) << -ry, rx, 1;

		Rot[i] = newRot*Rot[i]; //Safe for matrix-matrix multiplication with Eigen
		Rot[i] = compute_orthonormal(Rot[i]);
	}

}

Eigen::Matrix3d CurveDeformation::compute_orthonormal(Eigen::Matrix3d& rot) {
	Eigen::JacobiSVD<Eigen::MatrixXd> svd(rot, Eigen::ComputeFullU | Eigen::ComputeFullV);
	Eigen::Matrix3d VT = svd.matrixV().transpose();
	return svd.matrixU()*VT;
}

/*void CurveDeformation::final_L1_pos(Eigen::MatrixXd &V) {
	Eigen::VectorXd Bx(no_vertices + fixed_indices.size()), By(no_vertices + fixed_indices.size()), Bz(no_vertices + fixed_indices.size());
	Eigen::Vector3d Rl;
	for (int i = 0; i < no_vertices; i++) {
		Rl = Rot[i] * original_L1.row(i).transpose();
		Bx[i] = Rl(0);
		By[i] = Rl(1);
		Bz[i] = Rl(2);
	}

	for (int i = 0; i < fixed_indices.size(); i++) {
		Bx[no_vertices + i] = V(fixed_indices[i], 0) * CONSTRAINT_WEIGHT;
		By[no_vertices + i] = V(fixed_indices[i], 1) * CONSTRAINT_WEIGHT;
		Bz[no_vertices + i] = V(fixed_indices[i], 2) * CONSTRAINT_WEIGHT;
	}

	Eigen::VectorXd xpos = solverL1.solve(A_L1_T*Bx);
	Eigen::VectorXd ypos = solverL1.solve(A_L1_T*By);
	Eigen::VectorXd zpos = solverL1.solve(A_L1_T*Bz);
	for (int i = 0; i < no_vertices; i++) { //update the position of non-fixed vertices of the stroke that is being pulled on
		if (!is_fixed[i]) {
			V.row(vert_bindings[i]) << xpos[i], ypos[i], zpos[i];
		}
	}
	/*if(stroke_ID == 0) { //We're pulling on the original stroke
		for(int i = 0; i < no_vertices; i++) { //update the position of non-fixed vertices of the stroke that is being pulled on
			if(!is_fixed[i]) {
				V.row(vert_bindings[i]) << xpos[i], ypos[i], zpos[i];
			}
		}
	} else { //We're pulling on an added stroke, and we want to avoid deforming the original curve, so don't update positions of the original stroke
		for(int i = 0; i < no_vertices; i++) { //update the position of non-fixed stroke vertices
			if(!is_fixed[i] && !part_of_original_stroke[vert_bindings[i]]) {
				V.row(vert_bindings[i]) << xpos[i], ypos[i], zpos[i];
			}
		}
	}

}*/

void CurveDeformation::final_L1_pos_test(Eigen::MatrixXd &V, PulledCurve& curve) {
	Eigen::VectorXd Bx(curve.vertex_triplets.rows() + fixed_vertices_local.size()), By(curve.vertex_triplets.rows() + fixed_vertices_local.size()), Bz(curve.vertex_triplets.rows() + fixed_vertices_local.size());
	Eigen::Vector3d Rl;
	Eigen::Matrix3d r0, r1;
	for (int i = 0; i < curve.vertex_triplets.rows(); i++) {
		r0 = Rot[vertex_triplet_to_rot_idx[i][0]];
		r1 = Rot[vertex_triplet_to_rot_idx[i][1]];

		Rl = average_rot(r0, r1)* original_L1.row(i).transpose();
		Bx[i] = Rl(0);
		By[i] = Rl(1);
		Bz[i] = Rl(2);
	}

	for (int i = 0; i < curve.fixed_vertices.size(); i++) {
		Bx[curve.vertex_triplets.rows() + i] = V(curve.fixed_vertices[i], 0) *CONSTRAINT_WEIGHT;
		By[curve.vertex_triplets.rows() + i] = V(curve.fixed_vertices[i], 1) *CONSTRAINT_WEIGHT;
		Bz[curve.vertex_triplets.rows() + i] = V(curve.fixed_vertices[i], 2) *CONSTRAINT_WEIGHT;
	}

	Eigen::VectorXd xpos = solverL1.solve(A_L1_T*Bx);
	Eigen::VectorXd ypos = solverL1.solve(A_L1_T*By);
	Eigen::VectorXd zpos = solverL1.solve(A_L1_T*Bz);
	for (int i = 0; i < curve.vertices.size(); i++) { //update the position of non-fixed vertices of the stroke that is being pulled on
		if (!is_fixed[i]) {
			V.row(curve.vertices[i]) << xpos[i], ypos[i], zpos[i];
		}
	}
}

Eigen::Matrix3d CurveDeformation::average_rot(Eigen::Matrix3d& r0, Eigen::Matrix3d& r1) {
	Eigen::Matrix3d rnew;
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			rnew(i, j) = (r0(i, j) + r1(i, j)) / 2;
		}
		double length = rnew.row(i).norm();
		for (int j = 0; j < 3; j++) {
			rnew(i, j) /= length;
		}
	}
	return rnew;
}
