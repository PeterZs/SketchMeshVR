#include <stdlib.h>
#include <iostream>
#include <igl/readOFF.h>
#include <igl/viewer/Viewer.h>
#include <igl/vertex_triangle_adjacency.h>
#include <igl/adjacency_list.h>
#include <igl/per_face_normals.h>
#include <igl/per_vertex_normals.h>
#include <igl/per_corner_normals.h>
#include <igl/facet_components.h>
#include <igl/jet.h>
#include <igl/barycenter.h>
#include <cmath>
#include <igl/triangle_triangle_adjacency.h>
#include "main.h"

/*** insert any libigl headers here ***/

using namespace std;
//using System.Diagnostics;
using Viewer = igl::viewer::Viewer;

// Vertex array, #V x3
Eigen::MatrixXd V;
// Face array, #F x3
Eigen::MatrixXi F;
// Per-face normal array, #F x3
Eigen::MatrixXd FN;
// Per-vertex normal array, #V x3
Eigen::MatrixXd VN;
// Per-corner normal array, (3#F) x3
Eigen::MatrixXd CN;
// Vectors of indices for adjacency relations
std::vector<std::vector<int> > VF, VFi, VV;
// Integer vector of component IDs per face, #F x1
Eigen::VectorXi cid;
// Per-face color array, #F x3
Eigen::MatrixXd component_colors_per_face;

//Mouse interaction
enum ToolMode{DRAW, ADD, CUT, EXTRUDE, PULL, REMOVE, CHANGE, SMOOTH, NAVIGATE, NONE};
ToolMode tool_mode = NAVIGATE;
bool skip_standardcallback = false;
int down_mouse_x = -1, down_mouse_y = -1;

bool callback_key_down(Viewer& viewer, unsigned char key, int modifiers) {

	if (key == '1') {
		//viewer.data.clear();
		//viewer.data.set_mesh(V, F);
	}

	if (key == 'D') { //use capital letters
		//Draw initial curve/mesh
		tool_mode = DRAW;
	}
	if (key == 'P') {
		tool_mode = PULL;
	}
	if (key == 'N') {
		//Use navigation
		tool_mode = NAVIGATE;
	}

	//viewer.ngui->refresh(); //TODO: is this needed?
	return true;
}

bool callback_mouse_down(Viewer& viewer, int button, int modifier) {
	if (button == (int)Viewer::MouseButton::Right) {
		return false;
	}
	down_mouse_x = viewer.current_mouse_x;
	down_mouse_y = viewer.current_mouse_y;

	if (tool_mode == DRAW) { //Creating the first curve/mesh
		skip_standardcallback = true;
	}
	else if (tool_mode == PULL) { //Dragging an existing curve
		skip_standardcallback = true;
	}
	else if (tool_mode == NAVIGATE) { //Navigate through the screen
		skip_standardcallback = false; //We do want to use the navigation functionality
	}


	return skip_standardcallback; //Will make sure that we use standard navigation responses if we didn't do special actions and vice versa
}

bool callback_load_mesh(Viewer& viewer, string filename)
{
	igl::readOFF(filename, V, F);
	viewer.data.clear();
	viewer.data.set_mesh(V, F);
	viewer.data.compute_normals();
	viewer.core.align_camera_center(viewer.data.V);
	std::cout << filename.substr(filename.find_last_of("/") + 1) << endl;
	return true;
}

int main(int argc, char *argv[]) {
	// Show the mesh
	Viewer viewer;
	viewer.callback_key_down = callback_key_down;
	//viewer.callback_load_mesh = callback_load_mesh;
	

	if (argc == 2)
	{
		// Read mesh
		igl::readOFF(argv[1], V, F);

	}
	else
	{
		// Read mesh
		igl::readOFF("../data/plane.off", V, F);
	}
	callback_key_down(viewer, '1', 0);

	viewer.data.set_mesh(V, F);
	viewer.core.align_camera_center(V);
	viewer.launch();
}