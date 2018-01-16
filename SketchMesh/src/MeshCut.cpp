#include "MeshCut.h"
#include <iostream>
#include <igl/unproject.h>
#include <igl/edge_topology.h>
#include <igl/vertex_triangle_adjacency.h>
#include <igl/adjacency_list.h>
#include <igl/slice.h>
#include <igl/triangle/triangulate.h>
#include <string>
#include <sstream>
#include "LaplacianRemesh.h"

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
	cut_main(m, surface_path, stroke);
	
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
	Eigen::MatrixXi stroke_edges(boundary_vertices.rows(), 2);
	for (int i = 0; i < boundary_vertices.rows(); i++) {
		vec = m.V.row(boundary_vertices[i]) - center;
		boundary_vertices_2D.row(i) << vec.dot(x_vec), vec.dot(y_vec);
		stroke_edges.row(i) << i, ((i + 1) % boundary_vertices.size());
	}

	Eigen::MatrixXd V2;
	Eigen::MatrixXi F2;
	Eigen::MatrixXi vertex_markers, edge_markers;
	igl::triangle::triangulate(boundary_vertices_2D.leftCols(2), stroke_edges, Eigen::MatrixXd(0, 0), Eigen::MatrixXi::Constant(boundary_vertices_2D.rows(), 1, 1), Eigen::MatrixXi::Constant(stroke_edges.rows(), 1, 1), "Yq25", V2, F2, vertex_markers, edge_markers); //Capital Q silences triangle's output in cmd line. Also retrieves markers to indicate whether or not an edge/vertex is on the mesh boundary

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
			m.V.row(boundary_vertices[i]) = v_tmp.transpose();
		}
		else {
			Eigen::Vector3d v_tmp = center.transpose();
			v_tmp += x_vec*V2(i, 0);
			v_tmp += y_vec*V2(i, 1);
			m.V.row(original_v_size+i) << v_tmp.transpose(); //Add interior vertex of the cut plane to mesh
		}
	}

	int vert_idx_in_mesh, size_before = m.F.rows();
	m.F.conservativeResize(m.F.rows() + F2.rows(), Eigen::NoChange);
	for (int i = 0; i < F2.rows(); i++) {
		for (int j = 0; j < 3; j++) { //Go over the face vertices
			if (F2(i, j) < boundary_vertices.rows()) { //Original boundary vertex
				vert_idx_in_mesh = boundary_vertices[F2(i, j)];
			} else { //Interior to cut plane
				vert_idx_in_mesh = F2(i, j) + original_v_size; //Add index to the currently existing number of mesh vertices (which already includes the boundary vertices on the cut stroke)
			}
			m.F(size_before + i, 2 - j) = vert_idx_in_mesh; //To ensure right face orientation after implementing extrusion
		}
	}
	
	//TODO: see laplacianCut::cut lines 93-104 about sharp boundaries
}


void MeshCut::cut_main(Mesh& m, SurfacePath& surface_path, Stroke& stroke){
	Eigen::VectorXi boundary_vertices = LaplacianRemesh::remesh_cut_remove_inside(m, surface_path, stroke.viewer.core.model, stroke.viewer.core.view, stroke.viewer.core.proj, stroke.viewer.core.viewport);
	mesh_open_hole(boundary_vertices, m, stroke);
}


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
}

