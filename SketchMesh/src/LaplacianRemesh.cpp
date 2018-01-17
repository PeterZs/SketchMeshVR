#include "LaplacianRemesh.h"
#include <iostream>
#include <Eigen/Core>
#include <igl/edge_topology.h>
#include <igl/vertex_triangle_adjacency.h>
#include <igl/slice.h>
#include <igl/adjacency_list.h>
#include <igl/project.h>


using namespace std;
using namespace igl;

bool LaplacianRemesh::is_front_loop = true; //Determines extrusion or cut
bool LaplacianRemesh::remove_inside_faces = true;

Eigen::MatrixXi LaplacianRemesh::EV, LaplacianRemesh::FE, LaplacianRemesh::EF;
vector<vector<int>> LaplacianRemesh::VV;


Eigen::VectorXi LaplacianRemesh::remesh_cut_remove_inside(Mesh & m, SurfacePath & surface_path, Eigen::Matrix4f model, Eigen::Matrix4f view, Eigen::Matrix4f proj, Eigen::Vector4f viewport) {
	is_front_loop = false;
	remove_inside_faces = true;
	adjacency_list(m.F, VV);
	return remesh(m, surface_path, model, view, proj, viewport);
}

Eigen::VectorXi LaplacianRemesh::remesh_extrusion_remove_inside(Mesh & m, SurfacePath & surface_path, Eigen::Matrix4f model, Eigen::Matrix4f view, Eigen::Matrix4f proj, Eigen::Vector4f viewport) {
	is_front_loop = true;
	remove_inside_faces = true;
	adjacency_list(m.F, VV);
	return remesh(m, surface_path, model, view, proj, viewport);
}

Eigen::VectorXi LaplacianRemesh::remesh(Mesh& m, SurfacePath& surface_path, Eigen::Matrix4f model, Eigen::Matrix4f view, Eigen::Matrix4f proj, Eigen::Vector4f viewport) {
	Eigen::VectorXi result;
	vector<bool> dirty_face(m.F.rows());
	vector<int> dirty_vertices(m.V.rows());

	vector<vector<int>> VF, VI;
	igl::edge_topology(m.V, m.F, EV, FE, EF);
	igl::vertex_triangle_adjacency(m.V.rows(), m.F, VF, VI);

	int face = -1;
	int edge = -1;
	vector<PathElement> path = surface_path.get_path();
	if(!is_front_loop) { //Append first element at the end in case of cutting, to make it a loop and compatible with extrusion's loop
		path.push_back(path[0]);
	}

	//Find inside (+1) and outside (+1) vertices. TODO: check that this works and that +1 is inside and -1 is outside
	for(int i = 0; i < path.size(); i++) {
		if(path[i].get_type() == PathElement::FACE) {
			face = path[i].get_ID();
		} else {
			edge = path[i].get_ID();
			if(EF(edge, 0) == face) {
				dirty_vertices[EV(edge, 0)] += 1;
				dirty_vertices[EV(edge, 1)] -= 1;
			} else {
				dirty_vertices[EV(edge, 0)] -= 1;
				dirty_vertices[EV(edge, 1)] += 1;
			}
			face = (EF(edge, 0) == face) ? EF(edge, 1) : EF(edge, 0); //get the polygon on the other side of the edge
		}
	}

	//Collect vertices on the boundary
	vector<int> outer_boundary_vertices;
	vector<int> inner_boundary_vertices;
	for(int i = 0; i < m.V.rows(); i++) {
		if(dirty_vertices[i] < 0) {
			inner_boundary_vertices.push_back(i);
		} else if(dirty_vertices[i] > 0) {
			outer_boundary_vertices.push_back(i);
		}
	}

	//Collect faces along the path
	for(int i = 0; i < path.size(); i++) {
		if(path[i].get_type() == PathElement::EDGE) {
			dirty_face[EF(path[i].get_ID(), 0)] = true;
			dirty_face[EF(path[i].get_ID(), 1)] = true;
		}
	}

	if(remove_inside_faces) {
		//Collect all inside faces via recursive floodfill
		for(int i = 0; i < m.V.rows(); i++) {
			if(dirty_vertices[i] < 0) {
				for(int j = 0; j < VF[i].size(); j++) {
					if(!dirty_face[VF[i][j]]) {
						propagate_dirty_faces(VF[i][j], dirty_face);
					}
				}
			}
		}
	}

	//Remove dirty faces
	vector<int> clean_faces;
	for(int i = 0; i < m.F.rows(); i++) {
		if(!dirty_face[i]) {
			clean_faces.push_back(i);
		}
	}

	Eigen::MatrixXi tmp_F;
	Eigen::VectorXi row_idx, col_idx(3);
	row_idx = Eigen::VectorXi::Map(clean_faces.data(), clean_faces.size()); //Create an Eigen::VectorXi from a std::vector
	col_idx.col(0) << 0, 1, 2;
	igl::slice(m.F, row_idx, col_idx, tmp_F); //Keep only the clean faces in the mesh
	m.F = tmp_F;

	Eigen::MatrixXi EV_new;
	igl::edge_topology(m.V, m.F, EV, FE, EF); //TODO: NOTE: might need to save the EV that is outputted here in a new variable since maybe some methods might need the old EV
	igl::vertex_triangle_adjacency(m.V.rows(), m.F, VF, VI);


	//NOTE: Output from sort_boundary_vertices is not necessarily counter-clockwise
	outer_boundary_vertices = sort_boundary_vertices(path[0].get_vertex(), outer_boundary_vertices, m);
	if(!remove_inside_faces) {
		inner_boundary_vertices = sort_boundary_vertices(path[0].get_vertex(), inner_boundary_vertices, m);
	}

	//Do not use the last path vertex, as it is a copy of the first and creates unwanted behaviour
	vector<int> path_vertices;
	Eigen::MatrixX3d path_vert_positions(path.size() - 1, 3);
	int original_V_size = m.V.rows();
	for(int i = 0; i < path.size() - 1; i++) {
		path_vertices.push_back(original_V_size + i);
		path_vert_positions.row(i) << path[i].get_vertex().transpose();
	}

	//TODO: check in laplacianremesh line 650-656 if that's needed
	//TODO: check in laplacianremesh line 668-693 to see if we need to resample the stroke based on distance

	/*double unit_length = (is_front_loop) ? compute_average_distance_between_onPolygon_vertices(path, m) : compute_average_length_of_crossing_edges(path, m);
	if(path[0].get_vertex() != path[path.size() - 1].get_vertex()) { //Make a loop out of cut-strokes for similarity with extrusion method
		path_vert_positions.conservativeResize(path_vert_positions.rows() + 1, Eigen::NoChange);
		path_vert_positions.row(path_vert_positions.rows() - 1) << path[0].get_vertex().transpose();
	}

	if(get_length(path_vert_positions) / unit_length < 12) {
		unit_length = get_length(path_vert_positions) / 12.0;
	}

	Eigen::MatrixX3d tmp_path = resample_stroke(path_vert_positions); //Resamples the stroke by moving vertices to the middle of their neigbors --> implicitly results in smoothing/round stroke. not what we want
	*/
	int size_before = m.V.rows();
	m.V.conservativeResize(m.V.rows() + path.size() - 1, Eigen::NoChange);
	//TODO NOTE: NEED TO CHAGNE THIS TO USE THE OUTCOME OF RESAMPLE_BY_LENGTH...
	for(int i = 0; i < path.size() - 1; i++) {
		m.V.row(size_before + i) << path[i].get_vertex().transpose();
		//m.V.row(size_before + i) << tmp_path.row(i); //this version uses the smoothing but results in a roundish stroke.
	}



	Eigen::MatrixXd tmp_V;
	Eigen::VectorXi row_idx2, col_idx2(3);
	row_idx2 = Eigen::VectorXi::Map(outer_boundary_vertices.data(), outer_boundary_vertices.size()); //Create an Eigen::VectorXi from a std::vector
	col_idx2.col(0) << 0, 1, 2;
	igl::slice(m.V, row_idx2, col_idx2, tmp_V); //Keep only the boundary vertices in the mesh

	//Compute the mean of the boundary vertices on the side of the mesh that is removed
	//In the case of CUT this is used as a viewpoint to determine the orientation of the boundary vertices that remain
	Eigen::MatrixXd removed_V;
	Eigen::VectorXi row_idx_removed = Eigen::VectorXi::Map(inner_boundary_vertices.data(), inner_boundary_vertices.size());
	Eigen::VectorXi col_idx_removed(3);
	col_idx_removed.col(0) << 0, 1, 2;
	igl::slice(m.V, row_idx_removed, col_idx_removed, removed_V);
	Eigen::RowVector3d mean_viewpoint = removed_V.colwise().mean();


	Eigen::Matrix4f modelview = view * model;
	if(!is_counter_clockwise_boundaries(tmp_V, modelview, proj, viewport, mean_viewpoint, !is_front_loop)) {
		reverse(outer_boundary_vertices.begin(), outer_boundary_vertices.end());
	}


	stitch(path_vertices, outer_boundary_vertices, m);
	if(!remove_inside_faces) {
		reverse_path(path_vertices);
		stitch(path_vertices, inner_boundary_vertices, m);
	}

	//TODO: check if 708-728 in laplacianremesh is needed here

	return Eigen::VectorXi::Map(path_vertices.data(), path_vertices.size()); //Create an Eigen::VectorXi from a std::vector
}

void LaplacianRemesh::propagate_dirty_faces(int face, vector<bool>& dirty_face) {
	dirty_face[face] = true;
	for(int i = 0; i < 3; i++) {
		int next_face = (EF(FE(face, i), 0) == face) ? EF(FE(face, i), 1) : EF(FE(face, i), 0); //get an adjacent polygon
		if(!dirty_face[next_face]) {
			propagate_dirty_faces(next_face, dirty_face);
		}
	}
}

vector<int> LaplacianRemesh::sort_boundary_vertices(Eigen::Vector3d start_vertex, vector<int> boundary_vertices, Mesh& m) {
	vector<int> sorted_boundary_vertices;
	int start_v = find_closest(boundary_vertices, start_vertex, m);
	int v = start_v;
	int prev = start_v;

	Eigen::VectorXi tmp1, tmp2;
	int min_idx, max_idx, equal_pos, v_pos, max_val, only_option_idx;
	bool have_one_option = false;
	while(true) {
		sorted_boundary_vertices.push_back(v);
		for(int i = 0; i < VV[v].size(); i++) {
			if(std::find(boundary_vertices.begin(), boundary_vertices.end(), VV[v][i]) != boundary_vertices.end()) { //Only continue if the adjacent vertex is also a boundary vertex
				if(v < VV[v][i]) {
					min_idx = v;
					v_pos = 0;
					max_idx = VV[v][i];
				} else {
					min_idx = VV[v][i];
					v_pos = 1;
					max_idx = v;
				}
				tmp1 = EV.col(0).cwiseEqual(min_idx).cast<int>();
				tmp2 = EV.col(1).cwiseEqual(max_idx).cast<int>();
				max_val = (tmp1 + tmp2).maxCoeff(&equal_pos); //Tmp1 and tmp2 will contain a 1 on the positions where they equal the min and max_idx. When adding them, we'll get 2 on the row that contains both of them
				if((max_val != 2 && VV[v][i] != prev) || (max_val == 2 && EF(equal_pos, v_pos) == -1)) { //The first condition means that both vertices were only connected to faces that have been removed, so it is a valid move. The second condition means that the edge has a NULL face on the correct side of the edge, so proceed to adjacent vertex
					prev = v;
					v = VV[v][i];
					have_one_option = false; //Remove our backtracking option
					break;
				} else if(max_val != 2 && VV[v][i] == prev) { //Check if we have any other options before we start backtracking over a boundary
					only_option_idx = VV[v][i]; //Remember the vertex that we could backtrack to
					have_one_option = true;
					continue;
				}
			}
		}

		if(have_one_option) { //If we still have the backtracking vertex buffered, use it.
			prev = v;
			v = only_option_idx;
			have_one_option = false;
		}

		if(v == start_v) {
			break;
		}
	}
	return sorted_boundary_vertices;
}

/*vector<int> LaplacianRemesh::sort_boundary_vertices(Eigen::Vector3d start_vertex, std::vector<int> boundary_vertices, Mesh& m) {
	vector<int> sorted_boundary_vertices;
	int start_v = find_closest(boundary_vertices, start_vertex, m);
	int v = start_v;
	int prev = start_v;
	vector<int> visited(m.V.rows(), 0);
	bool found_another_one;
	while(true) {
		sorted_boundary_vertices.push_back(v);
		for(int i = 0; i < VV[v].size(); i++) {
			if((std::find(boundary_vertices.begin(), boundary_vertices.end(), VV[v][i]) != boundary_vertices.end()) && !visited[VV[v][i]]) { //Find a neighboring vertex that's also on the boundary, and hasn't been visited before (original boundary vertices are never set as visited, check in next line)
				if(VV[v][i] == prev) {
					if(m.part_of_original_stroke(VV[v][i])) {//The found vertex is on the original boundary and the previously visited vertex. First check if there's any more vertices available that are not on the original boundary, as the boundary vertex is the last point before wrapping around (has to do with loops)
						found_another_one = false;
						for(int j = i + 1; j < VV[v].size(); j++) {
							if((std::find(boundary_vertices.begin(), boundary_vertices.end(), VV[v][j]) != boundary_vertices.end()) && !visited[VV[v][j]] && VV[v][j] != prev) { //Another suitable vertex that ISN'T the previous one
								prev = v;
								v = VV[v][j];
								found_another_one = true;
								break;
							}
						}
						if(!found_another_one) {
							prev = v;
							v = VV[v][i];
						}

					} else { //The found vertex is the previously found vertex but isn't on the original boundary. Backtracking on interior vertices is not allowed, so search for next option!
						continue;
					}
				} else {
					prev = v;
					v = VV[v][i];
				}

				if(!m.part_of_original_stroke(prev) && prev != start_v) {
					visited[prev] = 1; //Do this afterwards, so that the start_v doesn't get set as visited (otherwise we will never reach it again and won't be able to break the loop)
				}
				break;
			}
		}
		if(v == start_v) {
			break;
		}
	}

	return sorted_boundary_vertices;
}
*/
int LaplacianRemesh::find_closest(vector<int> vertices, Eigen::Vector3d base, Mesh& m) {
	int closest = vertices[0];
	double min = (m.V.row(vertices[0]) - base.transpose()).norm();
	double d;

	for(int i = 1; i < vertices.size(); i++) {
		d = (m.V.row(vertices[i]) - base.transpose()).norm();
		if(d < min) {
			min = d;
			closest = vertices[i];
		}
	}
	return closest;
}

void LaplacianRemesh::stitch(std::vector<int> path_vertices, std::vector<int> boundary_vertices, Mesh& m) {
	int path_idx = 0;
	int outer_idx = 0;

	boundary_vertices = reorder(boundary_vertices, m.V.row(path_vertices[0]), m);

	Eigen::RowVector3d start_path_v = m.V.row(path_vertices[0]);
	Eigen::RowVector3d start_outer_v = m.V.row(boundary_vertices[0]);
	Eigen::RowVector3d path_v = start_path_v;
	Eigen::RowVector3d outer_v = start_outer_v;
	Eigen::RowVector3d next_path_v, next_outer_v;
	int path_v_idx = path_vertices[0], outer_v_idx = boundary_vertices[0], next_path_v_idx, next_outer_v_idx;

	bool proceed_outer_v;

	while(true) {
		next_path_v = m.V.row(path_vertices[(path_idx + 1) % path_vertices.size()]);
		next_outer_v = m.V.row(boundary_vertices[(outer_idx + 1) % boundary_vertices.size()]);
		next_path_v_idx = path_vertices[(path_idx + 1) % path_vertices.size()];
		next_outer_v_idx = boundary_vertices[(outer_idx + 1) % boundary_vertices.size()];

		proceed_outer_v = true;
		if(outer_idx == boundary_vertices.size()) {
			proceed_outer_v = false;
		} else if(path_idx == path_vertices.size()) {
			proceed_outer_v = true;
		} else if((path_v - next_outer_v).norm() < (outer_v - next_path_v).norm()) {
			proceed_outer_v = true;
		} else {
			proceed_outer_v = false;
		}

		//Add faces
		if(proceed_outer_v) {
			m.F.conservativeResize(m.F.rows() + 1, m.F.cols());
			m.F.row(m.F.rows() - 1) << path_v_idx, outer_v_idx, next_outer_v_idx;
			outer_v = next_outer_v;
			outer_v_idx = next_outer_v_idx;
			outer_idx++;
		} else {
			m.F.conservativeResize(m.F.rows() + 1, m.F.cols());
			m.F.row(m.F.rows() - 1) << next_path_v_idx, path_v_idx, outer_v_idx;
			path_v = next_path_v;
			path_v_idx = next_path_v_idx;
			path_idx++;
		}

		if(path_idx == path_vertices.size() && outer_idx == boundary_vertices.size()) {
			break;
		}
	}
}

vector<int> LaplacianRemesh::reorder(vector<int> boundary_vertices, Eigen::Vector3d start_v, Mesh& m) {
	int index = 0;
	double min = INFINITY;
	for(int i = 0; i < boundary_vertices.size(); i++) {
		Eigen::Vector3d vert = m.V.row(boundary_vertices[i]);
		double d = (vert - start_v).norm();
		if(d < min) {
			min = d;
			index = i;
		}
	}

	vector<int> reordered;
	for(int i = 0; i < boundary_vertices.size(); i++) {
		reordered.push_back(boundary_vertices[(i + index) % boundary_vertices.size()]);
	}
	return reordered;
}

//TODO: check this or remove if not necessary
void LaplacianRemesh::reverse_path(vector<int> path_vertices) {
	reverse(path_vertices.begin(), path_vertices.end());
	int v = path_vertices[path_vertices.size() - 1];
	path_vertices.insert(path_vertices.begin() + 1, path_vertices.begin(), path_vertices.end() - 1);
	path_vertices[0] = v;
}

double LaplacianRemesh::compute_average_length_of_crossing_edges(vector<PathElement> path, Mesh& m) {
	double total = 0.0;
	int count = 0;
	for(int i = 0; i < path.size(); i++) {
		if(path[i].get_type() == PathElement::EDGE) {
			total += (m.V.row(EV(path[i].get_ID(), 0)) - m.V.row(EV(path[i].get_ID(), 1))).norm();
			count++;
		}
	}
	return total / count;
}

double LaplacianRemesh::compute_average_distance_between_onPolygon_vertices(vector<PathElement> path, Mesh& m) {
	vector<int> face_vertices;
	for(int i = 0; i < path.size(); i++) {
		if(path[i].get_type() == PathElement::FACE) {
			face_vertices.push_back(i);
		}
	}

	double total = 0.0;
	Eigen::Vector3d v0, v1;
	for(int i = 0; i < face_vertices.size() - 1; i++) {
		v0 = path[face_vertices[i]].get_vertex();
		v1 = path[face_vertices[i + 1]].get_vertex();
		total += (v0 - v1).norm();
	}

	return total / (face_vertices.size() - 1);
}

double LaplacianRemesh::get_length(Eigen::MatrixX3d path_vertices) {
	double length = 0.0;
	Eigen::RowVector3d v0, v1;
	for(int i = 0; i < path_vertices.rows() - 1; i++) {
		v0 = path_vertices.row(i);
		v1 = path_vertices.row(i + 1);
		length += (v0 - v1).norm();
	}
	return length;
}

Eigen::MatrixX3d LaplacianRemesh::resample_stroke(Eigen::MatrixX3d & original_stroke3DPoints) {
	Eigen::MatrixX3d new_stroke3DPoints = Eigen::MatrixX3d::Zero(original_stroke3DPoints.rows(), 3);
	int nr_iterations = max(5.0, original_stroke3DPoints.rows() / 4.0);
	for(int i = 0; i < nr_iterations; i++) {
		move_to_middle(original_stroke3DPoints, new_stroke3DPoints);
		move_to_middle(new_stroke3DPoints, original_stroke3DPoints); //TODO: is this to save us a copy? basically performing an extra move_to_middle step
	}
	return new_stroke3DPoints;
}

void LaplacianRemesh::move_to_middle(Eigen::MatrixX3d &positions, Eigen::MatrixX3d &new_positions) {
	int n = positions.rows();
	Eigen::Vector3d prev, cur, next;

	for(int i = 0; i < n; i++) {
		prev = positions.row(((i - 1) + n) % n);
		cur = positions.row(i % n);
		next = positions.row(((i + 1) + n) % n);

		new_positions(i, 0) = (cur[0] * 2 + prev[0] + next[0]) / 4;
		new_positions(i, 1) = (cur[1] * 2 + prev[1] + next[1]) / 4;
		new_positions(i, 2) = (cur[2] * 2 + prev[2] + next[2]) / 4;
	}
}

bool LaplacianRemesh::is_counter_clockwise_boundaries(Eigen::MatrixXd boundary_points, Eigen::Matrix4f modelview, Eigen::Matrix4f proj, Eigen::Vector4f viewport, Eigen::RowVector3d mean_viewpoint, bool cut) {
	//First project the points to 2D
	Eigen::RowVector3d center = boundary_points.colwise().mean();
	Eigen::Vector3d normal(0, 0, 0);
	Eigen::Vector3d vec0, vec1;
	for(int i = 0; i < boundary_points.rows(); i++) {
		vec0 = boundary_points.row(i) - center;
		vec1 = boundary_points.row((i + 1) % boundary_points.rows()) - center;
		normal += vec1.cross(vec0);
	}
	normal.normalize();


	if(cut) {
		if(normal.dot(center - mean_viewpoint) < 0) {
			return false;
		} else {
			return true;
		}
	}

	double total_area = 0.0;
	boundary_points = boundary_points.rowwise() - center;
	Eigen::RowVector3d pt, vert;
	Eigen::Vector2d prev, next;
	vert = boundary_points.row(boundary_points.rows() - 1);
	igl::project(vert, modelview, proj, viewport, pt); //project the boundary vertex and store in pt
//	prev = vert.leftCols(2).transpose();
	prev = pt.leftCols(2).transpose();
	for(int i = 0; i < boundary_points.rows(); i++) {
		vert = boundary_points.row(i);
		igl::project(vert, modelview, proj, viewport, pt); //project the boundary vertex and store in pt
		//next = vert.leftCols(2).transpose();
		next = pt.leftCols(2).transpose();
		total_area += (prev[1] + next[1]) * (next[0] - prev[0]);
		prev = next;
	}

	if(total_area > 0) { //points are clockwise
		return false;
	}

	return true;
}