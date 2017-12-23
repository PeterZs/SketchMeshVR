#include "MeshCut.h"
#include <iostream>
#include <igl/unproject.h>
#include "Mesh.h"
#include "SurfacePath.h"

using namespace std;
using namespace igl;

int MeshCut::prev_vertex_count = -1;
int MeshCut::ID = -1;

void MeshCut::cut(Eigen::MatrixXd &V, Eigen::MatrixXi &F, Eigen::VectorXi &vertex_boundary_markers, Eigen::VectorXi &part_of_original_stroke, Stroke& stroke) {
	if(V.rows() != prev_vertex_count) {
		ID++;
		prev_vertex_count = V.rows();
	}
	Mesh m(V, F, vertex_boundary_markers, part_of_original_stroke, ID);

	SurfacePath surface_path;
	surface_path.create_from_stroke(stroke); //Prepares the drawn stroke (inserts extra points at the edges that it crosses)
	for(int i = 0; i < surface_path.get_path().size(); i++) { //NOTE: cut stroke points that lie inside the polygons aren't visible because they're "overwritten" by the pink ones
		cout << surface_path.get_path()[i].get_vertex().transpose() << endl;
		stroke.viewer.data.add_points(surface_path.get_path()[i].get_vertex().transpose(), Eigen::RowVector3d(1, 0, 0));
	}

	//cut_main(m);
}
