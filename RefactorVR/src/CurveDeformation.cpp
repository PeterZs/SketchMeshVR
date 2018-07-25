#include "CurveDeformation.h"
#include <igl/adjacency_list.h>
#include <igl/edge_topology.h>
#include <Eigen/Sparse>

using namespace std;
using namespace igl;

int moving_vertex_ID, prev_range_size = -1;
double prev_drag_size, DRAG_SCALE = 2.5;
Eigen::RowVector3d start_pos;
vector<vector<int>> CurveDeformation::neighbors;
Eigen::MatrixXi CurveDeformation::EV, CurveDeformation::FE, CurveDeformation::EF;

Eigen::VectorXi visited;
Eigen::VectorXd distance_to_vert;
Eigen::MatrixXd original_positions;
vector<CurveDeformation::PulledCurve> curves;

void CurveDeformation::startPullCurve(int _moving_vertex_ID, Eigen::MatrixXd& V, Eigen::MatrixXi& F) {
	moving_vertex_ID = _moving_vertex_ID;
	start_pos = V.row(moving_vertex_ID);
	prev_range_size = -1;
	prev_drag_size = 0.0;
	original_positions = V;
	adjacency_list(F, neighbors);
	igl::edge_topology(V, F, EV, FE, EF);
	visited.resize(V.rows());
	distance_to_vert.resize(V.rows());
}

void CurveDeformation::pullCurveTest(const Eigen::RowVector3d& pos, Eigen::MatrixXd& V, Eigen::VectorXi& edge_boundary_markers) {
	double drag_size = (pos - start_pos).norm() * DRAG_SCALE;
	bool ROI_is_updated = false;

	if (prev_drag_size < drag_size) { //Take the current drag_size and current_ROI_size relative to the size of the stroke we're pulling on. //TODO: test if we want to take it relative to the size of the whole mesh instead (we have access to V after all)?
		ROI_is_updated = update_ROI_test(drag_size, V, edge_boundary_markers);
	}
	V = original_positions; //TODO: see if we need this/whether it makes any difference in stability

	if (prev_range_size <= 1) {
		V.row(moving_vertex_ID) = pos;
	}
	else {

		if (ROI_is_updated) {
			for (int i = 0; i < curves.size(); i++) {
				curves[i].laplacian_curve_edit.setup_for_update_curve(curves[i].vertices, curves[i].fixed_vertices, curves[i].edges, curves[i].fixed_edges, curves[i].vertex_triplets, curves[i].edge_triplets, V, EV);
			}
		}
		V.row(moving_vertex_ID) = pos;
		for (int i = 0; i < curves.size(); i++) {
			curves[i].laplacian_curve_edit.update_curve(V);
		}

	}
}

bool CurveDeformation::update_ROI_test(double drag_size, Eigen::MatrixXd& V, Eigen::VectorXi& edge_boundary_markers) {
	prev_drag_size = drag_size;

	vector<int> vertices_in_range = collect_vertices_within_drag_length(drag_size, V, edge_boundary_markers);
	if (vertices_in_range.size() == prev_range_size) {
		return false;
	}
	prev_range_size = vertices_in_range.size();
	sort_by_distance(vertices_in_range);
	std::cout << "Vertices in range: ";
	for (int i = 0; i < vertices_in_range.size(); i++) {
		std::cout << vertices_in_range[i] << "   ";
	}
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
				if (fixed[curve.vertices[j]] && std::find(curve.fixed_vertices.begin(), curve.fixed_vertices.end(), curve.vertices[j]) == curve.fixed_vertices.end()) { //Fix vertices of a curve that are already part of a parent curve (hierarchy is based on vertex distance to handle vertex)
					curves.back().fixed_vertices.push_back(curve.vertices[j]);
				}
				fixed[curve.vertices[j]] = 1;
			}

			if (curves.back().fixed_vertices.size() == 1 && curves.back().vertices.size()>=3) { //Means that there is only 1 fixed vertex, which is the handle. All other vertices are free (so we'll end up dragging the mesh around instead of pulling on it)
				int count = 0;
				if (std::find(curves.back().vertices.begin(), curves.back().vertices.end(), vertices_in_range[vertices_in_range.size() - 1]) != curves.back().vertices.end() && std::find(curves.back().fixed_vertices.begin(), curves.back().fixed_vertices.end(), vertices_in_range[vertices_in_range.size()-1]) == curves.back().fixed_vertices.end()) { //If the last available vertex (furthest away) is part of this curve, but it is not fixed (means we have a loop?)
					curves.back().fixed_vertices.push_back(vertices_in_range[vertices_in_range.size() - 1]);
					count++;
				}
				if (std::find(curves.back().vertices.begin(), curves.back().vertices.end(), vertices_in_range[vertices_in_range.size() - 2]) != curves.back().vertices.end() && std::find(curves.back().fixed_vertices.begin(), curves.back().fixed_vertices.end(), vertices_in_range[vertices_in_range.size() - 2]) == curves.back().fixed_vertices.end()) {
					curves.back().fixed_vertices.push_back(vertices_in_range[vertices_in_range.size() - 2]);
					count++;
				}
				/*if (count == 2) {
					curves.back().fixed_edges.push_back(find_edge(curves.back().fixed_vertices[curves.back().fixed_vertices.size() - 1], curves.back().fixed_vertices[curves.back().fixed_vertices.size() - 2]));
				}*/
			}

			std::cout << "TEST" << std::endl;
			for (int i = 0; i < curves.back().fixed_vertices.size(); i++) {
				std::cout << curves.back().fixed_vertices[i] << std::endl;
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

	std::cout << "moving vertex id: " << moving_vertex_ID << std::endl;
	vector<int> vertices_in_range;
	vertices_in_range.push_back(moving_vertex_ID);
	propagate(-1, moving_vertex_ID, drag_size, vertices_in_range, V, edge_boundary_markers);
	
	return vertices_in_range;
}

void CurveDeformation::propagate(int prev_edge, int vert, double remaining_distance, std::vector<int>& vertices_in_range, Eigen::MatrixXd& V, Eigen::VectorXi& edge_boundary_markers) {
	int next_vert;
	double edge_length;
	for (int i = 0; i < neighbors[vert].size(); i++) {
		int edge = find_edge(vert, neighbors[vert][i]);
		if (edge == prev_edge) {
			continue;
		}
		edge_length = (original_positions.row(vert) - original_positions.row(neighbors[vert][i])).norm();
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

	int edge_marker = 0, edge;
	for (int i = 0; i < neighbors[vert].size(); i++) {
		edge = find_edge(vert, neighbors[vert][i]);
		if (edge_boundary_markers[edge] && !edge_consumed[edge] && visited[neighbors[vert][i]]) {
			edge_marker = edge_boundary_markers[edge];
		}
	}

	if (!edge_marker) {
		return false;
	}

	propagate(curve, -1, vert, edge_marker, edge_consumed, edge_boundary_markers);

	int prev_edge, next_edge;
	for (int i = 0; i < curve.edges.size(); i++) {
		prev_edge = get_adjacent_seam_edge(EV(curve.edges[i], 0), curve.edges[i], curve.edges, edge_boundary_markers);
		next_edge = get_adjacent_seam_edge(EV(curve.edges[i], 1), curve.edges[i], curve.edges, edge_boundary_markers);
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
	int edge;
	for (int i = 0; i < neighbors[vert].size(); i++) {
		edge = find_edge(vert, neighbors[vert][i]);
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
	int edge2;
	for (int i = 0; i < neighbors[vert].size(); i++) {
		edge2 = find_edge(vert, neighbors[vert][i]);
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
	int edge;
	for (int i = 0; i < neighbors[vert].size(); i++) {
		edge = find_edge(vert, neighbors[vert][i]);
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