#include "MeshExtrusion.h"
#include "LaplacianRemesh.h"
#include "Plane.h"
#include <iostream>
#include <algorithm>
#include <igl/unproject_ray.h>
#include <igl/triangle/triangulate.h>
#include <igl/edge_topology.h>
#include <igl/slice.h>

using namespace std;
using namespace igl;

bool MeshExtrusion::extrude_prepare(Stroke& base, SurfacePath& surface_path) {
	base.counter_clockwise();
	
	bool success = surface_path.create_from_stroke_extrude(base);

	if (!success) {
		return false;
	}
	post_extrude_prepare_update_points(base, surface_path);
	return true;
}

bool MeshExtrusion::extrude_main(Eigen::MatrixXd &V, Eigen::MatrixXi &F, Eigen::VectorXi &vertex_boundary_markers, Eigen::VectorXi &part_of_original_stroke, Eigen::VectorXi &new_mapped_indices, Eigen::VectorXi &sharp_edge, SurfacePath& surface_path, Stroke& stroke, Stroke& base, Eigen::Matrix4f model, Eigen::Matrix4f view, Eigen::Matrix4f proj, Eigen::Vector4f viewport) {

	Mesh m(V, F, vertex_boundary_markers, part_of_original_stroke, new_mapped_indices, sharp_edge, -1); //Give ID -1 since we don't use the mesh ID here anyway
	stroke.counter_clockwise();
	bool remesh_success = true;
	Eigen::VectorXi boundary_vertices = LaplacianRemesh::remesh_extrusion_remove_inside(m, surface_path, model, view, proj, viewport, remesh_success);

	if (!remesh_success) {
		//NOTE: remesh will have cleared m.new_mapped_indices but since this is cleared every time it is needed, there's no need to revert this
		return false;
	}

	vector<int> sharp_edge_indices;
	for(int i = 0; i < m.sharp_edge.rows(); i++) {
		if(m.sharp_edge[i]) {
			sharp_edge_indices.push_back(i);
		}
	}

	Eigen::MatrixXi EV, FE, EF;
	igl::edge_topology(m.V, m.F, EV, FE, EF);

    //Store the vertex indices for sharp edges in the original mesh
	Eigen::MatrixXi sharpEV;
	Eigen::VectorXi sharpEV_row_idx, sharpEV_col_idx(2);
	sharpEV_row_idx = Eigen::VectorXi::Map(sharp_edge_indices.data(), sharp_edge_indices.size());
	sharpEV_col_idx.col(0) << 0, 1;
	igl::slice(EV, sharpEV_row_idx, sharpEV_col_idx, sharpEV);


    //Compute the real 3D positions of the silhouette stroke
	Eigen::RowVector3d center(0, 0, 0);
	Eigen::Vector3d normal(0, 0, 0);
	get_normal_and_center(center, normal, m, boundary_vertices);

	Eigen::MatrixX3d sil_3D_points = stroke.get3DPoints();
	Eigen::MatrixXd silhouette_vertices(sil_3D_points.rows() - 1, 3);
	if (stroke.has_been_reversed) {
		silhouette_vertices = sil_3D_points.bottomRows(sil_3D_points.rows() - 1);
	} else {
		silhouette_vertices = sil_3D_points.topRows(sil_3D_points.rows()-1);
	}


	int most_left_vertex_idx = 0, most_right_vertex_idx = 0;
	find_left_and_right(most_left_vertex_idx, most_right_vertex_idx, m, boundary_vertices, silhouette_vertices.row(0), silhouette_vertices.row(silhouette_vertices.rows()-1), center);
	
	Eigen::RowVector3d dir = (silhouette_vertices.row(silhouette_vertices.rows()-1) - silhouette_vertices.row(0)).normalized();
	remove_out_of_bounds_silhouette(silhouette_vertices, center, m.V.row(boundary_vertices[most_left_vertex_idx]), m.V.row(boundary_vertices[most_right_vertex_idx]), dir);

	if(silhouette_vertices.rows() == 0) {
		cout << "No valid silhouette vertices. Try again." << endl;
		return false;
	}
	
	vector<int> sil_original_indices = add_silhouette_vertices(m, stroke.get_ID(), silhouette_vertices);
	vector<int> sil_original_indices_looped = sil_original_indices;
	sil_original_indices_looped.push_back(sil_original_indices[0]);
	stroke.set_closest_vert_bindings(sil_original_indices_looped);

	//Create the first loops
	Eigen::MatrixXd front_loop3D = silhouette_vertices;
	vector<int> front_loop_base_original_indices;
	create_loop(m, front_loop3D, boundary_vertices, front_loop_base_original_indices, most_left_vertex_idx, most_right_vertex_idx);
	//Create the second loop
	Eigen::MatrixXd back_loop3D = silhouette_vertices.colwise().reverse();
	vector<int> back_loop_base_original_indices;
	create_loop(m, back_loop3D, boundary_vertices, back_loop_base_original_indices, most_right_vertex_idx, most_left_vertex_idx);


	
	//Get center and normal for first half-dome and generate mesh
	Eigen::RowVector3d new_center_front(0,0,0);
	for (int i = 0; i < front_loop3D.rows(); i++) {
		new_center_front += front_loop3D.row(i);
	}
	for (int i = 0; i < silhouette_vertices.rows(); i++) {
		new_center_front -= silhouette_vertices.row(i);
	}
	new_center_front /= (front_loop3D.rows() - silhouette_vertices.rows());

	Eigen::Vector3d normal_front(0,0,0);
	Eigen::Vector3d vec0, vec1;
	for (int i = 0; i < front_loop3D.rows(); i++) {
		vec0 = front_loop3D.row(i) - new_center_front;
		vec1 = front_loop3D.row((i + 1) % front_loop3D.rows()) - new_center_front;
		normal_front += vec0.cross(vec1);
	}
	normal_front.normalize();

	Eigen::Vector3d x_vec = (m.V.row(boundary_vertices[most_right_vertex_idx]) - m.V.row(boundary_vertices[most_left_vertex_idx])).normalized();
	Eigen::Vector3d y_vec = normal_front.cross(x_vec);
	Eigen::Vector3d offset = x_vec.cross(y_vec);
	offset *= 0.05;

	generate_mesh(m, front_loop3D, new_center_front, x_vec, y_vec, offset, silhouette_vertices.rows(), front_loop_base_original_indices, sil_original_indices);

	//Get center and normal for second half-dome and generate mesh
	Eigen::RowVector3d new_center_back(0, 0, 0);
	for (int i = 0; i < back_loop3D.rows(); i++) {
		new_center_back += back_loop3D.row(i);
	}
	for (int i = 0; i < silhouette_vertices.rows(); i++) {
		new_center_back -= silhouette_vertices.row(i);
	}
	new_center_back /= (back_loop3D.rows() - silhouette_vertices.rows());


	Eigen::Vector3d normal_back(0, 0, 0);
	for (int i = 0; i < back_loop3D.rows(); i++) {
		vec0 = back_loop3D.row(i) - new_center_front;
		vec1 = back_loop3D.row((i + 1) % back_loop3D.rows()) - new_center_back;
		normal_back += vec0.cross(vec1);
	}
	normal_back.normalize();

	y_vec = normal_back.cross(x_vec);
	offset = x_vec.cross(y_vec);
	offset *= 0.05;

	reverse(sil_original_indices.begin(), sil_original_indices.end());
	generate_mesh(m, back_loop3D, new_center_back, x_vec, y_vec, offset, silhouette_vertices.rows(), back_loop_base_original_indices, sil_original_indices);

    //Update tracking variables with new vertex indices
	post_extrude_main_update_points(stroke, silhouette_vertices);
	post_extrude_main_update_bindings(base, surface_path);
	update_sharp_edges(m, sharpEV);
	return true;
}


void MeshExtrusion::remove_out_of_bounds_silhouette(Eigen::MatrixXd& silhouette_vertices, Eigen::RowVector3d& center, const Eigen::RowVector3d& left_most, const Eigen::RowVector3d& right_most, Eigen::RowVector3d& dir) {
	vector<int> keep_indices;
	Eigen::Vector3d center_to_vertex;
	double most_left_dot_prod = dir.dot(left_most - center);
	double most_right_dot_prod = dir.dot(right_most - center);
	double dot_prod;

	for (int i = 0; i < silhouette_vertices.rows(); i++) {
		center_to_vertex = silhouette_vertices.row(i) - center;
		dot_prod = dir.dot(center_to_vertex);
		if (dot_prod > most_left_dot_prod || dot_prod < most_right_dot_prod) { //Skip silhouette point because it's out of range
			continue;
		}
		else {
			keep_indices.push_back(i);
		}
	}

	Eigen::MatrixXd result;
	Eigen::VectorXi keep_row_idx, keep_col_idx(3);
	keep_row_idx = Eigen::VectorXi::Map(keep_indices.data(), keep_indices.size());
	keep_col_idx.col(0) << 0, 1, 2;
	igl::slice(silhouette_vertices, keep_row_idx, keep_col_idx, result);
	silhouette_vertices = result;
}

/** Generates a mesh inside the loop3D vertices. **/
void MeshExtrusion::generate_mesh(Mesh& m, Eigen::MatrixXd loop3D, Eigen::Vector3d center, Eigen::Vector3d x_vec, Eigen::Vector3d y_vec, Eigen::Vector3d offset, int nr_silhouette_vert, vector<int> loop_base_original_indices, vector<int> sil_original_indices) {
	int size_before_gen = m.V.rows();
	Eigen::MatrixXd loop2D(loop3D.rows(), 2);
	Eigen::Vector3d vec;
	Eigen::MatrixXi loop_stroke_edges(loop3D.rows(), 2);

	for(int i = 0; i < loop3D.rows(); i++) {
		vec = (loop3D.row(i).transpose() - center);
		loop2D.row(i) << vec.dot(x_vec), vec.dot(y_vec);
		loop_stroke_edges.row(i) << i, ((i + 1) % loop3D.rows());
	}

	Eigen::MatrixXd V2;
	Eigen::MatrixXi F2;
	Eigen::MatrixXi vertex_markers, edge_markers;
	igl::triangle::triangulate(loop2D, loop_stroke_edges, Eigen::MatrixXd(0, 0), Eigen::MatrixXi::Constant(loop2D.rows(), 1, 1), Eigen::MatrixXi::Constant(loop_stroke_edges.rows(), 1, 1), "Yq25", V2, F2, vertex_markers, edge_markers); //Capital Q silences triangle's output in cmd line. Also retrieves markers to indicate whether or not an edge/vertex is on the mesh boundary
							
	Eigen::RowVector3d vert;
	for(int i = 0; i < V2.rows(); i++) {
		if(i >= loop2D.rows()) { //Interior vertices
			vert = center;
			vert += x_vec.transpose()*V2(i, 0);
			vert += y_vec.transpose()*V2(i, 1);
			vert += offset.transpose();
			m.V.conservativeResize(m.V.rows() + 1, Eigen::NoChange);
			m.V.row(m.V.rows() - 1) = vert;
			m.part_of_original_stroke.conservativeResize(m.part_of_original_stroke.rows() + 1);
			m.vertex_boundary_markers.conservativeResize(m.vertex_boundary_markers.rows() + 1);
			m.part_of_original_stroke(m.part_of_original_stroke.rows() - 1) = 0;
			m.vertex_boundary_markers(m.vertex_boundary_markers.rows() - 1) = 0; //Interior mesh vertex, so non-boundary
		}
	}
	update_face_indices(m, F2, sil_original_indices, loop_base_original_indices, nr_silhouette_vert, size_before_gen, loop2D.rows());
	
}

/** Gets the center point in the extrusion base and the normal of its average plane. **/
void MeshExtrusion::get_normal_and_center(Eigen::RowVector3d& center, Eigen::Vector3d& normal, Mesh& m, Eigen::VectorXi &boundary_vertices) {
	for(int i = 0; i < boundary_vertices.rows(); i++) {
		center += m.V.row(boundary_vertices[i]);
	}
	center /= boundary_vertices.rows();

	Eigen::Vector3d vec0, vec1;
	for(int i = 0; i < boundary_vertices.rows(); i++) {
		vec0 = m.V.row(boundary_vertices[i]) - center;
		vec1 = m.V.row(boundary_vertices[(i + 1) % boundary_vertices.rows()]) - center;
		normal += vec0.cross(vec1);
	}
	normal.normalize();
}

void MeshExtrusion::find_left_and_right(int& most_left_vertex_idx, int& most_right_vertex_idx, Mesh& m, Eigen::VectorXi &boundary_vertices, const Eigen::RowVector3d& sil_start, const Eigen::RowVector3d& sil_end, Eigen::RowVector3d& center) {
	Eigen::Vector3d dir = sil_end - sil_start;
	dir.normalize();
	Eigen::Vector3d center_to_vertex;
	double dot_prod = dir.dot(m.V.row(boundary_vertices[0]) - center);
	double min = dot_prod, max = dot_prod;
	most_left_vertex_idx = 0;
	most_right_vertex_idx = 0;

	for (int i = 1; i < boundary_vertices.rows(); i++) {
		center_to_vertex = (m.V.row(boundary_vertices[i]) - center);
		dot_prod = dir.dot(center_to_vertex);
		if (dot_prod > max) {
			max = dot_prod;
			most_left_vertex_idx = i;
		}
		else if (dot_prod < min) {
			min = dot_prod;
			most_right_vertex_idx = i;
		}
	}
}

/** Adds the silhouette vertices to the mesh and also inserts their tracking variables. **/
vector<int> MeshExtrusion::add_silhouette_vertices(Mesh& m, int stroke_ID, Eigen::MatrixXd& silhouette_vertices) {
	vector<int> sil_original_indices;
	int size_before_silhouette = m.V.rows();
	m.V.conservativeResize(m.V.rows() + silhouette_vertices.rows(), Eigen::NoChange);
	m.part_of_original_stroke.conservativeResize(m.part_of_original_stroke.rows() + silhouette_vertices.rows());
	m.vertex_boundary_markers.conservativeResize(m.vertex_boundary_markers.rows() + silhouette_vertices.rows());
	for(int i = 0; i < silhouette_vertices.rows(); i++) {
		m.V.row(size_before_silhouette + i) = silhouette_vertices.row(i);
		m.part_of_original_stroke[size_before_silhouette + i] = 0;
		m.vertex_boundary_markers[size_before_silhouette + i] = stroke_ID;
		sil_original_indices.push_back(size_before_silhouette + i);
	}
	return sil_original_indices;
}

/** Creates a loop consisting of the silhouette vertices and half of the base stroke vertices. **/
void MeshExtrusion::create_loop(Mesh& m, Eigen::MatrixXd& loop3D, Eigen::VectorXi& boundary_vertices, vector<int> &loop_base_original_indices, int start_idx, int end_idx) {
	int idx = start_idx;
	while(true) {
		loop3D.conservativeResize(loop3D.rows() + 1, Eigen::NoChange);
		loop3D.row(loop3D.rows() - 1) = m.V.row(boundary_vertices[idx]);
		loop_base_original_indices.push_back(boundary_vertices[idx]);
		if(idx == end_idx) {
			break;
		}
		idx++;
		if(idx == boundary_vertices.rows()) {
			idx = 0;
		}
	}

}

/** Updates the sharp_edge tracker with new edge indices after the mesh topology changed. **/
void MeshExtrusion::update_sharp_edges(Mesh& m, Eigen::MatrixXi sharpEV) {
	Eigen::MatrixXi EV, FE, EF;
	igl::edge_topology(m.V, m.F, EV, FE, EF);
	m.sharp_edge.resize(EV.rows());
	m.sharp_edge.setZero();

	int start, end, equal_pos;
	Eigen::VectorXi col1Equals, col2Equals;
	for(int i = 0; i < sharpEV.rows(); i++) {
		start = sharpEV(i, 0); //No vertices were removed in the meantime, so use indices as-is
		end = sharpEV(i, 1);
		if(start == -1 || end == -1) { //Sharp edge no longer exists
			continue;
		}

		col1Equals = EV.col(0).cwiseEqual(std::min(start, end)).cast<int>();
		col2Equals = EV.col(1).cwiseEqual(std::max(start, end)).cast<int>();
		(col1Equals + col2Equals).maxCoeff(&equal_pos); //Find the row that contains both vertices of this edge

		m.sharp_edge[equal_pos] = 1; //Set this edge to be sharp
	}
}

/** Updates the old face indices with the new indexing and inserts new faces that are made in triangulation. **/
void MeshExtrusion::update_face_indices(Mesh& m, Eigen::MatrixXi& F2, vector<int> sil_original_indices, vector<int> loop_base_original_indices, int nr_silhouette_vert, int size_before_gen, int loop2D_size) {
	int original_F_size = m.F.rows();
	m.F.conservativeResize(m.F.rows() + F2.rows(), Eigen::NoChange);
	int vert_idx_in_mesh;
	for(int i = 0; i < F2.rows(); i++) {
		for(int j = 0; j < 3; j++) {
			if(F2(i, j) < nr_silhouette_vert) { //Silhouette vertex
				vert_idx_in_mesh = sil_original_indices[F2(i, j)];
			} else if(F2(i, j) < loop2D_size) { //Extrusion base vertex
				vert_idx_in_mesh = loop_base_original_indices[F2(i, j) - nr_silhouette_vert];
			} else {
				vert_idx_in_mesh = size_before_gen + F2(i, j) - loop2D_size; //Interior point
			}
			m.F(original_F_size + i, j) = vert_idx_in_mesh;
		}
	}
}

/** Updates the base stroke's 3DPoints with the new vertices. PathElements of the type FACE stay in the same position, but PathElements of the type EDGE get added. **/
void MeshExtrusion::post_extrude_prepare_update_points(Stroke& stroke, SurfacePath& surface_path) {
	vector<PathElement> path = surface_path.get_path();
	Eigen::MatrixX3d new_3DPoints(path.size(), 3);

	for(int i = 0; i < path.size(); i++) {
		new_3DPoints.row(i) = path[i].get_vertex().transpose();
	}

	stroke.set3DPoints(new_3DPoints);
}

/** Updates the silhouette stroke's 3DPoints with the new vertices. **/
void MeshExtrusion::post_extrude_main_update_points(Stroke &stroke, Eigen::MatrixXd new_positions) {
	Eigen::MatrixX3d new_3DPoints(new_positions.rows() + 1, 3);
	for(int i = 0; i < new_positions.rows(); i++) {
		new_3DPoints.row(i) = new_positions.row(i);
	}
	new_3DPoints.row(new_3DPoints.rows() - 1) = new_positions.row(0); //Make looped
	

	stroke.set3DPoints(new_3DPoints); 
}

/** Update the vertex bindings for the extrusion base stroke. Assumes indexing in m.V before any vertices are removed, so it will become the correct indexing after stroke.update_vert_bindings() is called. **/
void MeshExtrusion::post_extrude_main_update_bindings(Stroke& base, SurfacePath& surface_path) {
	vector<PathElement> path = surface_path.get_path();
	vector<int> new_closest_vertex_indices(path.size());
	for(int i = 0; i < path.size() - 1; i++) {
		new_closest_vertex_indices[i] = path[i].get_v_idx();
	}
	new_closest_vertex_indices[path.size() - 1] = new_closest_vertex_indices[0]; //In laplacianRemesh, the last path vertex is skipped because it is a duplicate of the first one. So need to set it manually

	base.set_closest_vert_bindings(new_closest_vertex_indices);
}
