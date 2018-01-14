#include "MeshExtrusion.h"
#include "LaplacianRemesh.h"
#include "Plane.h"
#include <iostream>
#include <igl/unproject_ray.h>
#include <igl/triangle/triangulate.h>

using namespace std;
using namespace igl;


int MeshExtrusion::prev_vertex_count = -1;
int MeshExtrusion::ID = -1;

void MeshExtrusion::extrude_prepare(Stroke& base, SurfacePath& surface_path) {
	base.counter_clockwise();

	//TODO: maybe need to do teddy.cleanstroke.clean like in Model line 1020 (ON BOTH BASE AND SILHOUETTE)
	surface_path.create_from_stroke_extrude(base);
}

void MeshExtrusion::extrude_main(Eigen::MatrixXd &V, Eigen::MatrixXi &F, Eigen::VectorXi &vertex_boundary_markers, Eigen::VectorXi &part_of_original_stroke, SurfacePath& surface_path, Stroke& stroke) {
	if(V.rows() != prev_vertex_count) {
		ID++;
		prev_vertex_count = V.rows();
	}

	Mesh m(V, F, vertex_boundary_markers, part_of_original_stroke, ID);
	stroke.counter_clockwise();
	Eigen::VectorXi boundary_vertices = LaplacianRemesh::remesh_extrusion_remove_inside(m, surface_path, stroke);

	//project points to 2D. TODO: for now following the example, should be able to work with libigl's project??
	Eigen::RowVector3d center(0, 0, 0);
	for(int i = 0; i < boundary_vertices.rows(); i++) {
		center += m.V.row(boundary_vertices[i]);
	}
	center /= boundary_vertices.rows();

	Eigen::Vector3d normal(0, 0, 0);
	Eigen::Vector3d vec0, vec1;
	for(int i = 0; i < boundary_vertices.rows(); i++) {
		vec0 = m.V.row(boundary_vertices[i]) - center;
		vec1 = m.V.row(boundary_vertices[(i + 1) % boundary_vertices.rows()]) - center;
		normal += vec1.cross(vec0);
	}
	normal.normalize();

	Eigen::Vector3d camera_to_center = center.transpose() - stroke.viewer.core.camera_eye.cast<double>();
	Eigen::Vector3d normal2 = normal.cross(camera_to_center);
	normal2.normalize();
	//Eigen::Vector3d normal2(1, 0, 0); //TODO: TEST THIS!! ABOVE IS WHAT TEDDY USES
	Eigen::Vector3d center_to_vertex = m.V.row(boundary_vertices[0]) - center;
	double dot_prod = normal2.dot(center_to_vertex);
	double max = dot_prod, min = dot_prod;
	int most_left_vertex_idx = 0, most_right_vertex_idx = 0;
	for(int i = 1; i < boundary_vertices.rows(); i++) {
		center_to_vertex = m.V.row(boundary_vertices[i]) - center;
		dot_prod = normal2.dot(center_to_vertex);
		if(dot_prod > max) {
			max = dot_prod;
			most_left_vertex_idx = i;
		} else if(dot_prod < min) {
			min = dot_prod;
			most_right_vertex_idx = i;
		}
	}

	Eigen::Vector3d v_tmp = (m.V.row(boundary_vertices[most_left_vertex_idx]) - m.V.row(boundary_vertices[most_right_vertex_idx])).transpose();
	Eigen::Vector3d norm_tmp = normal.cross(v_tmp);
	norm_tmp.normalize();
	Plane pop_surface(m.V.row(boundary_vertices[most_left_vertex_idx]), norm_tmp);

	Eigen::MatrixXd silhouette_vertices(0, 3);
	//TODO: CHECK WHAT IS GOING ON HERE. PROBABLY WON'T WORK LIKE THIS. MAYBE CREATE A PLANE THAT IS PERPENDICULAR TO THE BASE STROKE PLANE AND PROJECT THE SILHOUETTE STROKE ONTO THAT
	Eigen::RowVector3d v;
	Eigen::Vector3d source, dir;
	Eigen::Matrix4f modelview = stroke.viewer.core.view * stroke.viewer.core.model;
	for(int i = 0; i < stroke.get_stroke2DPoints().rows(); i++) {
		Eigen::Vector2d tmp = stroke.get_stroke2DPoints().row(i);
		unproject_ray(tmp, modelview, stroke.viewer.core.proj, stroke.viewer.core.viewport, source, dir);
		v = pop_surface.intersect_point(source, dir);
		center_to_vertex = v - center;
		dot_prod = normal2.dot(center_to_vertex);
		if(dot_prod > max || dot_prod < min) {
			cout << "silhouette point is outside of range " << i << endl;
			continue;
		}
		silhouette_vertices.conservativeResize(silhouette_vertices.rows() + 1, Eigen::NoChange);
		silhouette_vertices.row(silhouette_vertices.rows() - 1) = v;
	}

	//Possibly reverse silhouette_vertices, so they start closest to the most_left_vertex
	if((m.V.row(boundary_vertices[most_left_vertex_idx]) - silhouette_vertices.row(0)).norm() > (m.V.row(boundary_vertices[most_right_vertex_idx]) - silhouette_vertices.row(0)).norm()) {
		silhouette_vertices = silhouette_vertices.colwise().reverse().eval();
	}

	//TODO: FOR NOW SKIPPING OVER RESAMPLING IN LAPLACIANEXTRUSION LINE 311-315 (SEEMS TO BE A SIMPLE RESAMPLE THOUGH)
	int size_before_silhouette = m.V.rows();
	m.V.conservativeResize(m.V.rows() + silhouette_vertices.rows(), Eigen::NoChange);
	vector<int> sil_original_indices;
	for(int i = 0; i < silhouette_vertices.rows(); i++) {
		m.V.row(size_before_silhouette + i) = silhouette_vertices.row(i);
		sil_original_indices.push_back(size_before_silhouette + i);
	}

	//TODO: FOR NOW SKIPPING LINE 324-328 OF LAPLACIANEXTRUSION

	//Create one of the two loops
	Eigen::MatrixXd front_loop3D = silhouette_vertices;
	cout << "number of sil vertices: " << silhouette_vertices.rows() << " most_left_vertex_idx: " << most_left_vertex_idx << "  most_right_vertex_idx: "<< most_right_vertex_idx << endl;
	vector<int> front_loop_base_original_indices;
	int idx = most_right_vertex_idx;
	cout << "indices of vertices in front loop: ";
	while(true) {
		front_loop3D.conservativeResize(front_loop3D.rows() + 1, Eigen::NoChange);
		front_loop3D.row(front_loop3D.rows() - 1) = m.V.row(boundary_vertices[idx]);
		front_loop_base_original_indices.push_back(boundary_vertices[idx]);
		cout << boundary_vertices[idx] << " ";
		if(idx == most_left_vertex_idx) {
			break;
		}
		idx++;
		if(idx == boundary_vertices.rows()) {
			idx = 0;
		}
	}
	cout << endl;
	//Create the second loop
	Eigen::MatrixXd back_loop3D = silhouette_vertices.colwise().reverse();
	vector<int> back_loop_base_original_indices;
	int idx2 = most_left_vertex_idx;
	while(true) {
		back_loop3D.conservativeResize(back_loop3D.rows() + 1, Eigen::NoChange);
		back_loop3D.row(back_loop3D.rows() - 1) = m.V.row(boundary_vertices[idx2]);
		back_loop_base_original_indices.push_back(boundary_vertices[idx2]);
		if(idx2 == most_right_vertex_idx) {
			break;
		}
		idx2++;
		if(idx2 == boundary_vertices.rows()) {
			idx2 = 0;
		}
	}


	Eigen::Vector3d x_vec = m.V.row(boundary_vertices[most_right_vertex_idx]) - m.V.row(boundary_vertices[most_left_vertex_idx]);
	x_vec.normalize();
	Eigen::Vector3d y_vec = normal;
	Eigen::Vector3d offset = x_vec.cross(y_vec);
	offset *= 0.05;

	generate_mesh(m, front_loop3D, center, x_vec, y_vec, offset, silhouette_vertices.rows(), front_loop_base_original_indices, sil_original_indices);
	reverse(sil_original_indices.begin(), sil_original_indices.end());
	generate_mesh(m, back_loop3D, center, x_vec*-1, y_vec, offset*-1, silhouette_vertices.rows(), back_loop_base_original_indices, sil_original_indices);


	stroke.viewer.data.clear();
	stroke.viewer.data.set_mesh(m.V, m.F);
	Eigen::MatrixXd N_Faces;
	igl::per_face_normals(m.V, m.F, N_Faces);
	stroke.viewer.data.set_normals(N_Faces);


}

void MeshExtrusion::generate_mesh(Mesh& m, Eigen::MatrixXd loop3D, Eigen::Vector3d center, Eigen::Vector3d x_vec, Eigen::Vector3d y_vec, Eigen::Vector3d offset, int nr_silhouette_vert, vector<int> loop_base_original_indices, vector<int> sil_original_indices) {
	int size_before_gen = m.V.rows();
	//Generate mesh for loop
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
	cout << "3d input: " << loop3D << endl << endl;
	cout << "input to triangle: " << endl << loop2D << endl << endl;
	igl::triangle::triangulate(loop2D, loop_stroke_edges, Eigen::MatrixXd(0, 0), Eigen::MatrixXi::Constant(loop2D.rows(), 1, 1), Eigen::MatrixXi::Constant(loop_stroke_edges.rows(), 1, 1), "S", V2, F2, vertex_markers, edge_markers); //Capital Q silences triangle's output in cmd line. Also retrieves markers to indicate whether or not an edge/vertex is on the mesh boundary
																																																																	  
	Eigen::RowVector3d vert;
	for(int i = 0; i < V2.rows(); i++) {
		if(i < loop2D.rows()) {//Original boundary 
			vert = loop3D.row(i);
		} else {
			vert = center;
			vert += x_vec.transpose()*V2(i, 0);
			vert += y_vec.transpose()*V2(i, 1);
			vert += offset.transpose();
			m.V.conservativeResize(m.V.rows() + 1, Eigen::NoChange);
			m.V.row(m.V.rows() - 1) = vert;
		}
	}

	int original_F_size = m.F.rows();
	m.F.conservativeResize(m.F.rows() + F2.rows(), Eigen::NoChange);
	int vert_idx_in_mesh;
	for(int i = 0; i < F2.rows(); i++) {
		for(int j = 0; j < 3; j++) {
			if(F2(i, j) < nr_silhouette_vert) { //Silhouette vertex
				vert_idx_in_mesh = sil_original_indices[F2(i,j)];
			} else if(F2(i, j) < loop2D.rows()) { //Extrusion base vertex
				vert_idx_in_mesh = loop_base_original_indices[F2(i, j) - nr_silhouette_vert];
			} else {
				vert_idx_in_mesh = size_before_gen + F2(i, j) - loop2D.rows(); //Interior point
			}
			m.F(original_F_size + i, j) = vert_idx_in_mesh;
		}
	}
}