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
Eigen::VectorXi visited, distance_to_vert;
vector<CurveDeformation::PulledCurve> curves;
vector<vector<int>> vertex_triplet_to_rot_idx;
Eigen::VectorXi is_fixed;

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
	current_ROI_size = 0.0;
	adjacency_list(F, neighbors);
	igl::edge_topology(V, F, EV, FE, EF);
	visited.resize(V.rows());
	distance_to_vert.resize(V.rows());
}

//pos is the position to where the user dragged the vertex
/*void CurveDeformation::pullCurve(const Eigen::RowVector3d& pos, Eigen::MatrixXd& V) {
	double drag_size = (pos - start_pos).norm();
	drag_size /= curve_diag_length; //TODO: Not sure if this is neccessary? pUllmesh2 doesn't do it
	bool ROI_is_updated = false;
	if (!current_ROI_size || (fabs(drag_size - current_ROI_size) > 0.01)) { //Take the current drag_size and current_ROI_size relative to the size of the stroke we're pulling on. //TODO: test if we want to take it relative to the size of the whole mesh instead (we have access to V after all)?
		ROI_is_updated = update_ROI(drag_size);
	}

	if (no_ROI_vert == 0 || current_ROI_size == 0) { //If we have no other "free" vertices other than the handle vertex, we simply move it to the target position
		V.row(moving_vertex_ID) = pos;
	}
	else {
		if (ROI_is_updated) {
			setup_for_update_curve(V);
		}
		V.row(moving_vertex_ID) = pos;
		update_curve(V);
	}
}*/

double CurveDeformation::compute_curve_diag_length(Stroke& _stroke) {
	Eigen::Vector3d maxBB = _stroke.get3DPoints().colwise().maxCoeff();
	Eigen::Vector3d minBB = _stroke.get3DPoints().colwise().minCoeff();
	return (maxBB - minBB).norm();
}

/*bool CurveDeformation::update_ROI(double drag_size) {
	int no_ROI_vert_tmp;
	if (smooth_deform_mode && no_vertices > 1) {
		no_ROI_vert_tmp = max(round(0.17*no_vertices) + 1, min(round(drag_size * no_vertices) + 1, ceil(((no_vertices - 1) / 2) - 1))); //Determine how many vertices to the left and to the right to have free (at most half-1 of all vertices on each side, always at least 1/6th of the number of vertices + 1 vertex fixed, to take care of thin meshes with many vertices)
	}
	else {
		no_ROI_vert_tmp = max(0.0, min(round(drag_size / 4.0 * no_vertices), ceil(((no_vertices - 1) / 2) - 1))); //Determine how many vertices to the left and to the right to have free (at most half-1 of all vertices on each side)
	}

	if (((no_ROI_vert == no_ROI_vert_tmp) && prev_loop_type == stroke_is_loop) || no_ROI_vert_tmp == 0) { //number of vertices in ROI didn't change
		return false;
	}

	current_ROI_size = drag_size;
	no_ROI_vert = no_ROI_vert_tmp;
	prev_loop_type = stroke_is_loop;

	int ROI_1, ROI_2;
	compute_ROI_boundaries(ROI_1, ROI_2);
	setup_fixed_indices(ROI_1, ROI_2);
	return true;
}*/

void CurveDeformation::pullCurveTest(const Eigen::RowVector3d& pos, Eigen::MatrixXd& V, Eigen::VectorXi& edge_boundary_markers) {
	double drag_size = (pos - start_pos).norm();
	drag_size /= curve_diag_length; //TODO: Not sure if this is neccessary? pUllmesh2 doesn't do it
	bool ROI_is_updated = false;
	if (!current_ROI_size || (fabs(drag_size - current_ROI_size) > 0.01)) { //Take the current drag_size and current_ROI_size relative to the size of the stroke we're pulling on. //TODO: test if we want to take it relative to the size of the whole mesh instead (we have access to V after all)?
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
		return false;
	}
	prev_range_size = vertices_in_range.size();

	sort_by_distance(vertices_in_range);

	//Collect vertices in each each curve and construct cuve hierarchy
	Eigen::VectorXi edge_consumed(EV.rows()), fixed(V.rows());
	for (int i = 0; i < V.rows(); i++) {
		while (true) {
			PulledCurve curve;
			bool success = create_pulled_curve_by_propagation(i, edge_consumed, curve, edge_boundary_markers);
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
	visited[moving_vertex_ID] = 1;
	distance_to_vert[moving_vertex_ID] = drag_size;

	vector<int> vertices_in_range;
	vertices_in_range.push_back(moving_vertex_ID);
	propagate(-1, moving_vertex_ID, drag_size, vertices_in_range, V, edge_boundary_markers);

	return vertices_in_range;
}

void CurveDeformation::propagate(int prev_edge, int vert, double remaining_distance, std::vector<int>& vertices_in_range, Eigen::MatrixXd& V, Eigen::VectorXi& edge_boundary_markers) {
	int next_vert;
	for (int i = 0; i < neighbors[vert].size(); i++) {
		int edge = find_edge(vert, neighbors[vert][i]);
		if (edge == prev_edge) {
			continue;
		}
		double edge_length = (V.row(vert) - V.row(neighbors[vert][i])).norm();
		if (edge_boundary_markers[edge]) {
			if (remaining_distance - edge_length > 0) {
				next_vert = neighbors[vert][i];
				if (visited[next_vert]) {
					continue;
				}
				visited[next_vert] = true;
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

void CurveDeformation::sort_by_distance(std::vector<int> vertices_in_range) {
	sort(vertices_in_range.begin(), vertices_in_range.end(), distance_sorter);
}

bool CurveDeformation::distance_sorter(int a, int b) {
	return distance_to_vert[a] < distance_to_vert[b];
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

/*void CurveDeformation::compute_ROI_boundaries(int& ROI_1, int& ROI_2) {
	if (stroke_is_loop) {
		ROI_1 = (((handle_ID - no_ROI_vert) + no_vertices) % no_vertices);
		ROI_2 = (((handle_ID + no_ROI_vert) + no_vertices) % no_vertices);
	}
	else {
		ROI_1 = max(0, handle_ID - no_ROI_vert); //For non-loop strokes (e.g. added strokes), don't wrap around but instead cap at the first and last vertex of the stroke
		ROI_2 = min(no_vertices - 1, handle_ID + no_ROI_vert);
	}
}*/

/*void CurveDeformation::setup_fixed_indices(int ROI_1, int ROI_2) {
	vector<int> fixed, fixed_local;
	is_fixed = Eigen::VectorXi::Zero(no_vertices);

	if (ROI_1 < ROI_2) {
		for (int i = 0; i < ROI_1; i++) {
			is_fixed[i] = 1;
			fixed.push_back(vert_bindings[i]);
			fixed_local.push_back(i);
		}
		for (int i = ROI_2 + 1; i < no_vertices; i++) {
			is_fixed[i] = 1;
			fixed.push_back(vert_bindings[i]);
			fixed_local.push_back(i);
		}
	}
	else {
		for (int i = ROI_2 + 1; i < ROI_1; i++) {
			is_fixed[i] = 1;
			fixed.push_back(vert_bindings[i]);
			fixed_local.push_back(i);
		}
	}

	fixed.push_back(moving_vertex_ID);
	fixed_local.push_back(handle_ID);
	is_fixed[handle_ID] = 1;

	fixed_indices = Eigen::VectorXi::Map(fixed.data(), fixed.size());
	fixed_indices_local = Eigen::VectorXi::Map(fixed_local.data(), fixed_local.size());
}*/

/*void CurveDeformation::setup_for_update_curve(Eigen::MatrixXd& V) {
	A.resize(no_vertices * 3 + no_vertices * 9 + fixed_indices.size() * 3 + fixed_indices.size() * 3, no_vertices * 3 + no_vertices * 3); //Solve for x,y,z simultaneously since we cannot reuse A anyway
	B.resize(no_vertices * 3 + no_vertices * 9 + fixed_indices.size() * 3 + fixed_indices.size() * 3);
	original_L0.resize(no_vertices, 3);

	if (stroke_is_loop) {
		for (int i = 0; i < no_vertices; i++) {
			original_L0.row(i) = V.row(vert_bindings[i]) - V.row(((vert_bindings[((i - 1) + no_vertices) % no_vertices]))); //This assumes that the stroke is looped, which might not always be true for added control strokes.
		}
	}
	else {
		original_L0.row(0) = V.row(vert_bindings[0]) - V.row(vert_bindings[1]);
		for (int i = 1; i < no_vertices; i++) {
			original_L0.row(i) = V.row(vert_bindings[i]) - V.row(vert_bindings[i - 1]);
		}

	}
	setup_for_L1_position_step(V);
}*/

void CurveDeformation::setup_for_update_curve_test(PulledCurve& curve, Eigen::MatrixXd& V) {
	A.resize(curve.edges.size() * 3 + curve.edge_triplets.rows() * 9 + curve.fixed_vertices.size() * 3 + curve.fixed_edges.size() * 3, curve.vertices.size() * 3 + curve.edges.size() * 3);
	B.resize(curve.edges.size() * 3 + curve.edge_triplets.rows() * 9 + curve.fixed_vertices.size() * 3 + curve.fixed_edges.size() * 3);
	Rot.resize(curve.edges.size());
	original_L0.resize(curve.edges.size(), 3);
	for (int i = 0; i < curve.edges.size(); i++) {
		original_L0.row(i) = V.row(EV(curve.edges[i],0)) - V.row(EV(curve.edges[i],1));
	}

	is_fixed.resize(curve.vertices.size());
	for (int i = 0; i < curve.fixed_vertices.size(); i++) {
		is_fixed[i] = 1;
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
		vertex_triplet_to_rot_idx[i][0] = find_edge(v0, v1);
		vertex_triplet_to_rot_idx[i][1] = find_edge(v1, v2);

		tripletList.push_back(T(i, v0, -0.5));
		tripletList.push_back(T(i, v0, 1.0));
		tripletList.push_back(T(i, v0, -0.5));
	}

	for (int i = 0; i < curve.fixed_vertices.size(); i++) {
		tripletList.push_back(T(curve.vertex_triplets.rows() + i, curve.fixed_vertices[i], CONSTRAINT_WEIGHT));
	}
	Eigen::SparseMatrix<double> A_L1(curve.vertex_triplets.rows() + curve.fixed_vertices.size(), curve.vertices.size());
	A_L1.setFromTriplets(tripletList.begin(), tripletList.end());

	A_L1_T = A_L1.transpose();
	solverL1.compute(A_L1_T*A_L1);
}

/*void CurveDeformation::setup_for_L1_position_step(Eigen::MatrixXd& V) {
	Eigen::SparseMatrix<double> A_L1(no_vertices + fixed_indices.size(), no_vertices);
	std::vector<T> tripletList;
	tripletList.reserve(no_vertices * 3);

	original_L1.resize(no_vertices, 3);
	int cur, next, prev;
	if (stroke_is_loop) {
		for (int i = 0; i < no_vertices; i++) {
			cur = vert_bindings[i];
			next = vert_bindings[(i + 1) % no_vertices];
			prev = vert_bindings[(((i - 1) + no_vertices) % no_vertices)];

			tripletList.push_back(T(i, (((i - 1) + no_vertices) % no_vertices), -0.5));
			tripletList.push_back(T(i, i, 1));
			tripletList.push_back(T(i, (i + 1) % no_vertices, -0.5));

			original_L1.row(i) = V.row(cur) - (0.5*V.row(next) + 0.5*V.row(prev));
		}
	}
	else {
		tripletList.push_back(T(0, 0, 1));
		tripletList.push_back(T(0, 1, -1));
		original_L1.row(0) = V.row(vert_bindings[0]) - V.row(vert_bindings[1]);

		tripletList.push_back(T(no_vertices - 1, no_vertices - 1, 1));
		tripletList.push_back(T(no_vertices - 1, no_vertices - 2, -1));

		original_L1.row(no_vertices - 1) = V.row(vert_bindings[no_vertices - 1]) - V.row(vert_bindings[no_vertices - 2]);

		for (int i = 1; i < no_vertices - 1; i++) {
			cur = vert_bindings[i];
			next = vert_bindings[(i + 1) % no_vertices];
			prev = vert_bindings[(((i - 1) + no_vertices) % no_vertices)];

			tripletList.push_back(T(i, i - 1, -0.5));
			tripletList.push_back(T(i, i, 1));
			tripletList.push_back(T(i, i + 1, -0.5));

			original_L1.row(i) = V.row(cur) - (0.5*V.row(next) + 0.5*V.row(prev));
		}
	}

	for (int i = 0; i < fixed_indices.size(); i++) {
		tripletList.push_back(T(no_vertices + i, fixed_indices_local[i], CONSTRAINT_WEIGHT));
	}

	A_L1.setFromTriplets(tripletList.begin(), tripletList.end());

	A_L1_T = A_L1.transpose();
	solverL1.compute(A_L1_T*A_L1);
}*/

/*void CurveDeformation::update_curve(Eigen::MatrixXd& V) {
	for (int i = 0; i < no_vertices; i++) {
		Rot[i] = Eigen::Matrix3d::Identity();
	}

	for (int i = 0; i < 2; i++) {
		solve_for_pos_and_rot(V);
		update_rot();
	}
	final_L1_pos(V);
}*/

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
	tripletList.reserve(curve.edges.size()*15 + curve.edge_triplets.rows()*3*27 + curve.fixed_vertices.size()*3 + curve.fixed_edges.size() * 3);

	int v0, v1;
	for (int i = 0; i < curve.edges.size(); i++) { //Laplacian of the position: v1' - v0' = dR * R * L(v1 - v0)
		v0 = EV(curve.edges[i], 0);
		v1 = EV(curve.edges[i], 1);
		
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
		e0 = curve.edge_triplets(i, 0);
		e1 = curve.edge_triplets(i, 1); 
		e2 = curve.edge_triplets(i, 2);

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

	for (int i = 0; i < curve.fixed_vertices.size(); i++) { //Position of fixed vertices: v_i - v_i' = 0
		tripletList.push_back(T(curve.edges.size() * 3 + curve.edge_triplets.rows() * 9 + i * 3, curve.fixed_vertices[i] * 3, CONSTRAINT_WEIGHT));
		tripletList.push_back(T(curve.edges.size() * 3 + curve.edge_triplets.rows() * 9 + i * 3 + 1, curve.fixed_vertices[i] * 3 + 1, CONSTRAINT_WEIGHT));
		tripletList.push_back(T(curve.edges.size() * 3 + curve.edge_triplets.rows() * 9 + i * 3 + 2, curve.fixed_vertices[i] * 3 + 2, CONSTRAINT_WEIGHT));

		B[curve.edges.size() * 3 + curve.edge_triplets.rows() * 9 + i * 3] = CONSTRAINT_WEIGHT * V(curve.fixed_vertices[i], 0);
		B[curve.edges.size() * 3 + curve.edge_triplets.rows() * 9 + i * 3 + 1] = CONSTRAINT_WEIGHT * V(curve.fixed_vertices[i], 1);
		B[curve.edges.size() * 3 + curve.edge_triplets.rows() * 9 + i * 3 + 2] = CONSTRAINT_WEIGHT * V(curve.fixed_vertices[i], 2);
	}

	for (int i = 0; i < curve.fixed_edges.size(); i++) {
		tripletList.push_back(T(curve.edges.size() * 3 + curve.edge_triplets.rows() * 9 + curve.fixed_vertices.size() * 3 + i * 3 + 0, curve.vertices.size() * 3 + curve.fixed_edges[i] * 3, CONSTRAINT_WEIGHT));
		tripletList.push_back(T(curve.edges.size() * 3 + curve.edge_triplets.rows() * 9 + curve.fixed_vertices.size() * 3 + i * 3 + 1, curve.vertices.size() * 3 + curve.fixed_edges[i] * 3 + 1, CONSTRAINT_WEIGHT));
		tripletList.push_back(T(curve.edges.size() * 3 + curve.edge_triplets.rows() * 9 + curve.fixed_vertices.size() * 3 + i * 3 + 2, curve.vertices.size() * 3 + curve.fixed_edges[i] * 3 + 2, CONSTRAINT_WEIGHT));

		B[curve.edges.size() * 3 + curve.edge_triplets.rows() * 9 + curve.fixed_vertices.size() * 3 + i * 3 + 0, curve.vertices.size() * 3 + curve.fixed_edges[i] * 3] = 0;
		B[curve.edges.size() * 3 + curve.edge_triplets.rows() * 9 + curve.fixed_vertices.size() * 3 + i * 3 + 0, curve.vertices.size() * 3 + curve.fixed_edges[i] * 3 + 1] = 0;
		B[curve.edges.size() * 3 + curve.edge_triplets.rows() * 9 + curve.fixed_vertices.size() * 3 + i * 3 + 0, curve.vertices.size() * 3 + curve.fixed_edges[i] * 3 + 2] = 0;
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
	Eigen::VectorXd Bx(curve.vertex_triplets.rows() + curve.fixed_vertices.size()), By(curve.vertex_triplets.rows() + curve.fixed_vertices.size()), Bz(curve.vertex_triplets.rows() + curve.fixed_vertices.size());
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
