#include "MeshExtrusion.h"
#include "LaplacianRemesh.h"
#include <iostream>

using namespace std;
using namespace igl;


int MeshExtrusion::prev_vertex_count = -1;
int MeshExtrusion::ID = -1;

void MeshExtrusion::extrude(Eigen::MatrixXd &V, Eigen::MatrixXi &F, Eigen::VectorXi &vertex_boundary_markers, Eigen::VectorXi &part_of_original_stroke, Stroke& base, Stroke& silhouette) {
	if(V.rows() != prev_vertex_count) {
		ID++;
		prev_vertex_count = V.rows();
	}

	Mesh m(V, F, vertex_boundary_markers, part_of_original_stroke, ID);
	base.counter_clockwise();
	silhouette.counter_clockwise();

	//TODO: maybe need to do teddy.cleanstroke.clean like in Model line 1020 (ON BOTH BASE AND SILHOUETTE)
	SurfacePath surface_path;
	surface_path.create_from_stroke_extrude(base);
	base.viewer.data.add_points(surface_path.get_path()[1].get_vertex().transpose(), Eigen::RowVector3d(0, 0, 0));

	extrude_main(m, surface_path, silhouette);

	base.viewer.data.clear();
	base.viewer.data.set_mesh(m.V, m.F);
	Eigen::MatrixXd N_Faces;
	igl::per_face_normals(m.V, m.F, N_Faces);
	base.viewer.data.set_normals(N_Faces);
}

void MeshExtrusion::extrude_main(Mesh& m, SurfacePath& surface_path, Stroke& stroke) {
	Eigen::VectorXi boundary_vertices = LaplacianRemesh::remesh_extrusion_remove_inside(m, surface_path);

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


	Eigen::Vector3d camera_to_center = center - stroke.viewer.core.camera_eye.cast<double>();
	Eigen::Vector3d normal2 = normal.cross(camera_to_center);

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
		}else if(dot_prod < min) {
			min = dot_prod;
			most_right_vertex_idx = i;
		}
	}

	cout << "test " << stroke.viewer.core.camera_center << endl << endl << stroke.viewer.core.camera_eye << endl;
	//Eigen::Vector3d cam_to_center = center - stroke.viewer.core.camera_eye;
	cout << m.V << endl;

}