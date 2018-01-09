#include "MeshCut.h"
#include <iostream>
#include <igl/unproject.h>
#include <igl/edge_topology.h>
#include <igl/vertex_triangle_adjacency.h>
#include <igl/adjacency_list.h>
#include <igl/slice.h>
#include <igl/triangle/triangulate.h>
#include "Mesh.h"
#include "SurfacePath.h"
#include <string>
#include <sstream>

using namespace std;
using namespace igl;

int MeshCut::prev_vertex_count = -1;
int MeshCut::ID = -1;
Eigen::MatrixXi MeshCut::EV, MeshCut::FE, MeshCut::EF;
vector<vector<int>> VV;

void MeshCut::cut(Eigen::MatrixXd &V, Eigen::MatrixXi &F, Eigen::VectorXi &vertex_boundary_markers, Eigen::VectorXi &part_of_original_stroke, Stroke& stroke) {
	if(V.rows() != prev_vertex_count) {
		ID++;
		prev_vertex_count = V.rows();
	}
	Mesh m(V, F, vertex_boundary_markers, part_of_original_stroke, ID);
	SurfacePath surface_path;
	surface_path.create_from_stroke(stroke); //Prepares the drawn stroke (inserts extra points at the edges that it crosses)
	adjacency_list(m.F, VV);
	cout << "So goo2" << endl;
	cut_main(m, surface_path, stroke);
	cout << "doing stuff" << endl;
	
	stroke.viewer.data.clear();
	stroke.viewer.data.set_mesh(m.V, m.F);
	Eigen::MatrixXd N_Faces;
	igl::per_face_normals(m.V, m.F, N_Faces);
	stroke.viewer.data.set_normals(N_Faces);
	return;
}


void MeshCut::mesh_open_hole(Eigen::VectorXi& boundary_vertices, Mesh& m, Stroke& stroke) {
	//project points to 2D. TODO: for now following the example, should be able to work with libigl's project??
	Eigen::RowVector3d center(0,0,0);
	for (int i = 0; i < boundary_vertices.rows(); i++) {
		center += m.V.row(boundary_vertices[i]);
	}
	center /= boundary_vertices.rows();

	Eigen::Vector3d normal(0, 0, 0);
	Eigen::Vector3d vec0, vec1;
	for (int i = 0; i < boundary_vertices.rows(); i++) {
		vec0 = m.V.row(boundary_vertices[i]) - center;
		vec1 = m.V.row(boundary_vertices[(i + 1) % boundary_vertices.rows()]) - center;
		normal += vec1.cross(vec0);
	}
	normal.normalize();

	Eigen::Vector3d x_vec = m.V.row(boundary_vertices[0]) - center;
	x_vec.normalize();
	Eigen::Vector3d y_vec = normal.cross(x_vec);

	Eigen::MatrixXd boundary_vertices_2D(boundary_vertices.rows(), 2);
	Eigen::Vector3d vec;
	Eigen::MatrixXi stroke_edges;
	stroke_edges.resize(boundary_vertices.rows(), 2);
	for (int i = 0; i < boundary_vertices.rows(); i++) {
		vec = m.V.row(boundary_vertices[i]) - center;
		boundary_vertices_2D.row(i) << vec.dot(x_vec), vec.dot(y_vec);
		stroke_edges.row(i) << i, ((i + 1) % boundary_vertices.size());
	}

	Eigen::Matrix4f modelview = stroke.viewer.core.view * stroke.viewer.core.model;
	Eigen::MatrixXd points_to_project;
	igl::slice(m.V, boundary_vertices, 1, points_to_project);
	//igl::project(points_to_project, modelview, stroke.viewer.core.proj, stroke.viewer.core.viewport, boundary_vertices_2D);
	Eigen::MatrixXd V2;
	Eigen::MatrixXi F2;
	Eigen::MatrixXi vertex_markers, edge_markers;
	igl::triangle::triangulate(boundary_vertices_2D.leftCols(2), stroke_edges, Eigen::MatrixXd(0, 0), Eigen::MatrixXi::Constant(boundary_vertices_2D.rows(), 1, 1), Eigen::MatrixXi::Constant(stroke_edges.rows(), 1, 1), "q25", V2, F2, vertex_markers, edge_markers); //Capital Q silences triangle's output in cmd line. Also retrieves markers to indicate whether or not an edge/vertex is on the mesh boundary

	int original_v_size = m.V.rows()-boundary_vertices.rows();		
	m.V.conservativeResize(original_v_size + V2.rows(), Eigen::NoChange);
	//project back to 3D
	//TODO: Additional Steiner points may be added between boundary vertices. These HAVE to be unprojected to 3D because we don't have an original 3D position. However, if we don't do the same with the original boundary vertices, we will get holes in the mesh. Applying the same unprojection to the original
	//boundary vertices results in loss of the drawn stroke shape
	for (int i = 0; i < V2.rows(); i++) {
		if (i < boundary_vertices.rows()) {//Original boundary 
			Eigen::Vector3d v_tmp = center.transpose();
			v_tmp += x_vec*V2(i, 0);
			v_tmp += y_vec*V2(i, 1);
		//	m.V.row(original_v_size+i) = v_tmp.transpose();
			m.V.row(boundary_vertices[i]) = v_tmp.transpose();
		}


		else {
			Eigen::Vector3d v_tmp = center.transpose();
			v_tmp += x_vec*V2(i, 0);
			v_tmp += y_vec*V2(i, 1);
			m.V.row(original_v_size+i) << v_tmp.transpose(); //Add interior vertex of the cut plane to mesh
		}
	}

	int vert_idx_in_mesh;
	for (int i = 0; i < F2.rows(); i++) {
		m.F.conservativeResize(m.F.rows() + 1, Eigen::NoChange);
		for (int j = 0; j < 3; j++) { //Go over the face vertices
			if (F2(i, j) < boundary_vertices.rows()) { //Original boundary vertex
				vert_idx_in_mesh = boundary_vertices[F2(i, j)];
			} else { //Interior to cut plane
				vert_idx_in_mesh = F2(i, j) + original_v_size; //Add index to the currently existing number of mesh vertices (which already includes the boundary vertices on the cut stroke)
			}
			m.F(m.F.rows() - 1, j) = vert_idx_in_mesh;
		}
	}
	
	//TODO: see laplacianCut::cut lines 93-104 about sharp boundaries
}


void MeshCut::cut_main(Mesh& m, SurfacePath& surface_path, Stroke& stroke){
	Eigen::VectorXi boundary_vertices = remesh_cut_remove_inside(m, surface_path);
	cout << "Done with 1" << endl;
	mesh_open_hole(boundary_vertices, m, stroke);
}

Eigen::VectorXi MeshCut::remesh_cut_remove_inside(Mesh& m, SurfacePath& surface_path) {
	bool is_front_loop = false;
	bool remove_inside_faces = true;
	Eigen::VectorXi result;

	vector<bool> dirty_face(m.F.rows());
	vector<int> dirty_vertices(m.V.rows());

	vector<vector<int>> VF, VI;
	igl::edge_topology(m.V, m.F, EV, FE, EF);
	igl::vertex_triangle_adjacency(m.V.rows(), m.F, VF, VI);


	int face = -1;
	int edge = -1;
	vector<PathElement> path = surface_path.get_path();

	//Find inside (+1) and outside (+1) vertices. TODO: check that this works and that +1 is inside and -1 is outside
	for (int i = 0; i < path.size(); i++) {
		if (path[i].get_type() == PathElement::FACE) {
			face = path[i].get_ID();
		} else {
			edge = path[i].get_ID();
			if (EF(edge,0)==face) {
				dirty_vertices[EV(edge, 0)] += 1;
				dirty_vertices[EV(edge, 1)] -= 1;
			}
			else {
				dirty_vertices[EV(edge, 0)] -= 1;
				dirty_vertices[EV(edge, 1)] += 1;
			}
			face = (EF(edge, 0) == face) ? EF(edge, 1) : EF(edge, 0); //get the polygon on the other side of the edge
		}
	}

	//Collect vertices on the boundary
	vector<int> outer_boundary_vertices;
	vector<int> inner_boundary_vertices;
	for (int i = 0; i < m.V.rows(); i++) {
		if (dirty_vertices[i] > 0) {
			inner_boundary_vertices.push_back(i);
		}
		else if (dirty_vertices[i] < 0) {
			outer_boundary_vertices.push_back(i);
		}
	}


	//Collect faces along the path
	for (int i = 0; i < path.size(); i++) {
		if (path[i].get_type() == PathElement::EDGE) {
			dirty_face[EF(path[i].get_ID(), 0)] = true;
			dirty_face[EF(path[i].get_ID(), 1)] = true;
		}
	}

	if (remove_inside_faces) {
		//Collect all inside faces via recursive floodfill
		for (int i = 0; i < m.V.rows(); i++) {
			if (dirty_vertices[i] > 0) {
				for (int j = 0; j < VF[i].size(); j++) {
					if (!dirty_face[VF[i][j]]) {
						propagate_dirty_faces(VF[i][j], dirty_face);
					}
				}
			}
		}
	}

	//Remove dirty faces
	vector<int> clean_faces;
	for (int i = 0; i < m.F.rows(); i++) {
		if (!dirty_face[i]) {
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
	igl::edge_topology(m.V, m.F, EV_new, FE, EF);
	igl::vertex_triangle_adjacency(m.V.rows(), m.F, VF, VI);
	cout << "getting out" << endl;

	outer_boundary_vertices = sort_boundary_vertices(surface_path.get_path()[0].get_vertex(), outer_boundary_vertices, m);
	if (!remove_inside_faces) {
		inner_boundary_vertices = sort_boundary_vertices(surface_path.get_path()[0].get_vertex(), inner_boundary_vertices, m);
	}
	cout << "getting out" << endl;

	vector<int> path_vertices;
	Eigen::MatrixX3d path_vert_positions(path.size(), 3);
	int original_V_size = m.V.rows();
	for (int i = 0; i < path.size(); i++) {
		path_vertices.push_back(original_V_size+i);
		path_vert_positions.row(i) << path[i].get_vertex().transpose();
	}

	//TODO: check in laplacianremesh line 650-656 if that's needed


	//TODO: check in laplacianremesh line 668-693 to see if we need to resample the stroke based on distance


	double unit_length = (is_front_loop) ? compute_average_distance_between_onPolygon_vertices(path, m) : compute_average_length_of_crossing_edges(path, m);
	if (path[0].get_vertex() != path[path.size()-1].get_vertex()) { //Make a loop out of cut-strokes for similarity with extrusion method
		path_vert_positions.conservativeResize(path_vert_positions.rows() + 1, Eigen::NoChange);
		path_vert_positions.row(path_vert_positions.rows() - 1) << path[0].get_vertex().transpose(); 
	}

	if(get_length(path_vert_positions) / unit_length < 12) {
		unit_length = get_length(path_vert_positions) / 12.0;
	}

	Eigen::MatrixX3d tmp_path = resample_stroke(path_vert_positions); //Resamples the stroke by moving vertices to the middle of their neigbors --> implicitly results in smoothing/round stroke. not what we want

	int size_before = m.V.rows();
	m.V.conservativeResize(m.V.rows() + path.size(), Eigen::NoChange);
	//TODO NOTE: NEED TO CHAGNE THIS TO USE THE OUTCOME OF RESAMPLE_BY_LENGTH...
	for(int i = 0; i < path.size(); i++) {
		m.V.row(size_before + i) << path[i].get_vertex().transpose();
		//m.V.row(size_before + i) << tmp_path.row(i); //this version uses the smoothing but results in a roundish stroke.
	}

	stitch(path_vertices, outer_boundary_vertices, m);
	if (!remove_inside_faces) {
		reverse_path(path_vertices);
		stitch(path_vertices, inner_boundary_vertices, m);
	}


	//TODO: check if 708-728 in laplacianremesh is needed here

	return Eigen::VectorXi::Map(path_vertices.data(), path_vertices.size()); //Create an Eigen::VectorXi from a std::vector
}


void MeshCut::propagate_dirty_faces(int face, vector<bool>& dirty_face) {
	dirty_face[face] = true;
	for (int i = 0; i < 3; i++) {
		int next_face = (EF(FE(face, i), 0) == face) ? EF(FE(face, i), 1) : EF(FE(face, i), 0); //get an adjacent polygon
		if (!dirty_face[next_face]) {
			propagate_dirty_faces(next_face, dirty_face);
		}
	}
}

vector<int> MeshCut::sort_boundary_vertices(Eigen::Vector3d start_vertex, std::vector<int> boundary_vertices, Mesh& m) {
	vector<int> sorted_boundary_vertices;
	int start_v = find_closest(boundary_vertices, start_vertex, m);
	int v = start_v;
	int prev = start_v;
	vector<int> visited(m.V.rows(), 0);
	cout <<  m.part_of_original_stroke << endl;
	bool found_another_one;
	while (true) {
		sorted_boundary_vertices.push_back(v);
		for (int i = 0; i < VV[v].size(); i++) {
			if ((std::find(boundary_vertices.begin(), boundary_vertices.end(), VV[v][i]) != boundary_vertices.end()) && !visited[VV[v][i]]) { //Find a neighboring vertex that's also on the boundary, and hasn't been visited before (original boundary vertices are never set as visited, check in next line)
				if(VV[v][i] == prev) {
					if(m.part_of_original_stroke(VV[v][i])) {//The found vertex is on the original boundary and the previously visited vertex. First check if there's any more vertices available that are not on the original boundary, as the boundary vertex is the last point before wrapping around (has to do with loops)
						//do stuff
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

				if(!m.part_of_original_stroke(prev)) {
					visited[prev] = 1; //Do this afterwards, so that the start_v doesn't get set as visited (otherwise we will never reach it again and won't be able to break the loop)
				}
				break;
			}
		}
		if (v == start_v) {
			break;
		}
	}
	cout << " yippie" << endl;
	return sorted_boundary_vertices;
}

 int MeshCut::find_closest(vector<int> vertices, Eigen::Vector3d base, Mesh& m) {
	 int closest = vertices[0];
	 double min = (m.V.row(vertices[0]) - base.transpose()).norm(); 
	 double d;

	 for (int i = 1; i<vertices.size(); i++) {
		d = (m.V.row(vertices[i]) - base.transpose()).norm(); 
		if (d < min) {
			min = d;
			closest = vertices[i];
		}
	}
	return closest;
}

void MeshCut::stitch(std::vector<int> path_vertices, std::vector<int> boundary_vertices, Mesh& m) {
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

	while (true) {
		next_path_v = m.V.row(path_vertices[(path_idx + 1) % path_vertices.size()]);
		next_outer_v = m.V.row(boundary_vertices[(outer_idx + 1) % boundary_vertices.size()]);
		next_path_v_idx = path_vertices[(path_idx + 1) % path_vertices.size()];
		next_outer_v_idx = boundary_vertices[(outer_idx + 1) % boundary_vertices.size()];
		
		proceed_outer_v = true;
		if (outer_idx == boundary_vertices.size()) {
			proceed_outer_v = false;
		}
		else if (path_idx == path_vertices.size()) {
			proceed_outer_v = true;
		}
		else if ((path_v - next_outer_v).norm() < (outer_v - next_path_v).norm()) { 
			proceed_outer_v = true;
		}
		else {
			proceed_outer_v = false;
		}

		//Add faces
		if (proceed_outer_v) {
			m.F.conservativeResize(m.F.rows() + 1, m.F.cols());
			m.F.row(m.F.rows() - 1) << path_v_idx, next_outer_v_idx, outer_v_idx;
			outer_v = next_outer_v;
			outer_v_idx = next_outer_v_idx;
			outer_idx++;
		}
		else {
			m.F.conservativeResize(m.F.rows() + 1, m.F.cols());
			m.F.row(m.F.rows() - 1) << next_path_v_idx, outer_v_idx, path_v_idx;
			path_v = next_path_v;
			path_v_idx = next_path_v_idx;
			path_idx++;
		}

		if (path_idx == path_vertices.size() && outer_idx == boundary_vertices.size()) {
			break;
		} 
	}
}

vector<int> MeshCut::reorder(vector<int> boundary_vertices, Eigen::Vector3d start_v, Mesh& m) {
	int index = 0;
	double min = INFINITY;
	for (int i = 0; i < boundary_vertices.size(); i++) {
		Eigen::Vector3d vert = m.V.row(boundary_vertices[i]);
		double d = (vert-start_v).norm();
		if (d < min) {
			min = d;
			index = i;
		}
	}

	vector<int> reordered;
	for (int i = 0; i < boundary_vertices.size(); i++) {
		reordered.push_back(boundary_vertices[(i + index) % boundary_vertices.size()]);
	}
	return reordered;
}

/*void MeshCut::reverse_path(Eigen::MatrixXd path_vertices) {
	path_vertices.colwise().reverse();
	Eigen::RowVector3d v = path_vertices.row(path_vertices.rows() - 1);
	path_vertices.block(1, 0, path_vertices.rows() - 1, path_vertices.cols()) = path_vertices.block(0, 0, path_vertices.rows() - 1, path_vertices.cols());
	path_vertices.row(0) = v;
}*/

void MeshCut::reverse_path(vector<int> path_vertices) {
	reverse(path_vertices.begin(), path_vertices.end());
	int v = path_vertices[path_vertices.size()];
	path_vertices.insert(path_vertices.begin() + 1, path_vertices.begin(), path_vertices.end() - 1);
	path_vertices[0] = v;
}

double MeshCut::compute_average_length_of_crossing_edges(vector<PathElement> path, Mesh& m) {
	double total = 0.0;
	int count = 0;
	for(int i = 0; i < path.size(); i++) {
		if(path[i].get_type() == PathElement::EDGE) {
			cout << EV.row(path[i].get_ID()) << endl;
			total += (m.V.row(EV(path[i].get_ID(), 0)) - m.V.row(EV(path[i].get_ID(), 1))).norm();
			count++;
		}
	}
	return total / count;
}

double MeshCut::compute_average_distance_between_onPolygon_vertices(vector<PathElement> path, Mesh& m) {
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
		v1 = path[face_vertices[i+1]].get_vertex();
		total += (v0 - v1).norm();
	}

	return total / (face_vertices.size() - 1);
}


double MeshCut::get_length(Eigen::MatrixX3d path_vertices) {
	double length = 0.0;
	Eigen::RowVector3d v0, v1;
	for(int i = 0; i < path_vertices.rows() - 1; i++) {
		v0 = path_vertices.row(i);
		v1 = path_vertices.row(i + 1);
		length += (v0 - v1).norm();
	}
	return length;
}
/*
Eigen::MatrixXd MeshCut::resample_by_length_with_fixes(vector<int> path_vertices, Mesh& m, double unit_length) {
	if(path_vertices.size() <= 1) {
		return m.V.row(path_vertices[0]);
	}

	//TODO: Skipping stuff from CleanStroke3D line 64-76. IS IT NEEDED?

	Eigen::MatrixXd resampled(0, 3);
	resampled.conservativeResize(resampled.rows() + 1, Eigen::NoChange);
	resampled.row(resampled.rows() - 1) = m.V.row(path_vertices[0]);

	int idx0 = 0, idx1;
	while(true) {
	//	idx1 = find_next_fi
	}
}*/



Eigen::MatrixX3d MeshCut::resample_stroke(Eigen::MatrixX3d & original_stroke3DPoints) {
	Eigen::MatrixX3d new_stroke3DPoints = Eigen::MatrixX3d::Zero(original_stroke3DPoints.rows(), 3);
	int nr_iterations = max(5.0, original_stroke3DPoints.rows() / 4.0);
	for(int i = 0; i < nr_iterations; i++) {
		move_to_middle(original_stroke3DPoints, new_stroke3DPoints);
		move_to_middle(new_stroke3DPoints, original_stroke3DPoints); //TODO: is this to save us a copy? basically performing an extra move_to_middle step
	}
	return new_stroke3DPoints;
}

void MeshCut::move_to_middle(Eigen::MatrixX3d &positions, Eigen::MatrixX3d &new_positions) {
	int n = positions.rows();
	//Do seperately for i=0, because modulo gives -1
	Eigen::Vector3d prev = positions.row(n - 1);
	Eigen::Vector3d cur = positions.row(0);
	Eigen::Vector3d next = positions.row(1);
	new_positions(0, 0) = (cur[0] * 2 + prev[0] + next[0]) / 4;
	new_positions(0, 1) = (cur[1] * 2 + prev[1] + next[1]) / 4;
	new_positions(0, 2) = (cur[2] * 2 + prev[2] + next[2]) / 4;

	for(int i = 1; i < n; i++) {
		prev = positions.row((i - 1) % n);
		cur = positions.row(i%n);
		next = positions.row((i + 1) % n);

		new_positions(i, 0) = (cur[0] * 2 + prev[0] + next[0]) / 4;
		new_positions(i, 1) = (cur[1] * 2 + prev[1] + next[1]) / 4;
		new_positions(i, 2) = (cur[2] * 2 + prev[2] + next[2]) / 4;
	}
}