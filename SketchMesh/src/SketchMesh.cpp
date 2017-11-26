#ifdef _WIN32
#include <Windows.h>
#else
#include <unistd.h>
#endif
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
#include <igl/cat.h>
#include <cmath>
#include <igl/triangle_triangle_adjacency.h>
#include "SketchMesh.h"
#include "Stroke.h"
#include "SurfaceSmoothing.h"
#include "CurveDeformation.h"

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

// Per vertex indicator of whether vertex is on boundary (on boundary if == 1)
Eigen::VectorXi vertex_boundary_markers;

//Mouse interaction
enum ToolMode{DRAW, ADD, CUT, EXTRUDE, PULL, REMOVE, CHANGE, SMOOTH, NAVIGATE, NONE};
ToolMode tool_mode = NAVIGATE;
bool skip_standardcallback = false;
int down_mouse_x = -1, down_mouse_y = -1;

double vertex_weights;
int smooth_iter = 1;
bool mouse_is_down = false; //We need this due to mouse_down not working in the nanogui menu, whilst mouse_up does work there

//For selecting vertices
//std::unique_ptr<Stroke> _stroke;
Stroke* _stroke;
int handleID = -1;

bool callback_key_down(Viewer& viewer, unsigned char key, int modifiers) {
	if (key == '1') {
		viewer.data.clear();
	}else if (key == 'D') { //use capital letters
		//Draw initial curve/mesh
		tool_mode = DRAW;
		_stroke->strokeReset();
	}else if (key == 'P') {
		tool_mode = PULL;
	}else if (key == 'N') {
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
	mouse_is_down = true;

	down_mouse_x = viewer.current_mouse_x;
	down_mouse_y = viewer.current_mouse_y;

	if (tool_mode == DRAW) { //Creating the first curve/mesh
		viewer.data.clear();
		_stroke->strokeReset();
		_stroke->strokeAddSegment(down_mouse_x, down_mouse_y);
		skip_standardcallback = true;
	}
	else if (tool_mode == PULL) { //Dragging an existing curve
		handleID = _stroke->selectClosestVertex(down_mouse_x, down_mouse_y);
		CurveDeformation::startPullCurve(*_stroke, handleID);
		skip_standardcallback = true;
	}
	else if (tool_mode == NAVIGATE) { //Navigate through the screen
		skip_standardcallback = false; //We do want to use the navigation functionality
	}
	else if(tool_mode == EXTRUDE) {
		_stroke->strokeAddSegmentExtrusion(down_mouse_x, down_mouse_y);
		skip_standardcallback = true;
	}


	return skip_standardcallback; //Will make sure that we use standard navigation responses if we didn't do special actions and vice versa
}


bool callback_mouse_move(Viewer& viewer, int mouse_x, int mouse_y) {
	if (!skip_standardcallback) {
		return false;
	}
	if (tool_mode == DRAW && viewer.down) { //If we're still holding the mouse down
		_stroke->strokeAddSegment(mouse_x, mouse_y);
		return true;
	} else if(tool_mode == EXTRUDE && viewer.down) {
		_stroke->strokeAddSegmentExtrusion(mouse_x, mouse_y);
		return true;
	} else if(tool_mode == PULL && viewer.down) {
		double x = mouse_x;
		double y = viewer.core.viewport(3) - mouse_y;

		Eigen::Matrix4f modelview = viewer.core.view * viewer.core.model;
		Eigen::RowVector3f pt1(viewer.data.V(handleID, 0), viewer.data.V(handleID, 1), viewer.data.V(handleID, 2));
		Eigen::RowVector3f pr;
		igl::project(pt1, modelview, viewer.core.proj, viewer.core.viewport, pr);
		Eigen::RowVector3d pt = igl::unproject(Eigen::Vector3f(x, y, pr[2]), modelview, viewer.core.proj, viewer.core.viewport).transpose().cast<double>();


		CurveDeformation::pullCurve(pt, V);
		viewer.data.set_mesh(V, F);
		viewer.data.compute_normals();
	}
	return false;
}

bool callback_mouse_up(Viewer& viewer, int button, int modifier) {
	if(!mouse_is_down) {
		return true;
	}
	mouse_is_down = false;

	if(tool_mode == DRAW) {
		if(_stroke->toLoop()) {//Returns false if the stroke only consists of 1 point (user just clicked)
            //Give some time to show the stroke
            #ifdef _WIN32
                Sleep(200);
            #else
                usleep(200000);  /* sleep for 200 milliSeconds */
            #endif
            _stroke->generate3DMeshFromStroke(vertex_boundary_markers);
			F = viewer.data.F;
			V = viewer.data.V;

			for(int i = 0; i < smooth_iter; i++) {
				SurfaceSmoothing::smooth(V, F, vertex_boundary_markers);
			}

			viewer.data.set_mesh(V, F);
            viewer.data.compute_normals();

            //Overlay the drawn stroke
			int strokeSize = (vertex_boundary_markers.array() > 0).count();
			Eigen::MatrixXd strokePoints = V.block(0, 0, strokeSize, 3);
			Eigen::MatrixXd tmp_1 = V.block(1, 0, strokeSize - 1, 3);
			Eigen::MatrixXd tmp_2 = V.row(0);
			Eigen::MatrixXd endPoints = igl::cat(1, tmp_1, tmp_2);
			viewer.data.add_points(strokePoints, Eigen::RowVector3d(1, 0, 0));
			viewer.data.add_edges(V.block(0, 0, strokeSize, 3), endPoints, Eigen::RowVector3d(1,0,0));
		}
		skip_standardcallback = false;
	}
	return skip_standardcallback;
}

//TODO: make callback for this in viewer, like in exercise 5 of shapemod
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
	viewer.callback_mouse_down = callback_mouse_down;
	viewer.callback_mouse_move = callback_mouse_move;
	viewer.callback_mouse_up = callback_mouse_up;
	//viewer.callback_load_mesh = callback_load_mesh;
    
    viewer.callback_init = [&](igl::viewer::Viewer& viewer)
    {
        // Add new group
        viewer.ngui->addGroup("Inflation");
        
        // Expose a variable directly ...
        viewer.ngui->addVariable("Vertex Weights",SurfaceSmoothing::vertex_weight);
        viewer.ngui->addVariable("Edge Weights",SurfaceSmoothing::edge_weight);

        
        // Expose a variable directly ...
        viewer.ngui->addVariable("Smoothing iterations",smooth_iter);

        
        // Add a button
        viewer.ngui->addButton("Perform 1 smoothing iteration",[&viewer](){
            SurfaceSmoothing::smooth(V,F,vertex_boundary_markers);
            viewer.data.set_mesh(V, F);
            viewer.data.compute_normals();
        });
        
        // call to generate menu
        viewer.screen->performLayout();
        return false;
    };

	//Init stroke selector
	//_stroke = std::unique_ptr<Stroke>(new Stroke(V, F, viewer));
	_stroke = new Stroke(V, F, viewer);
	if (argc == 2)
	{
		// Read mesh
		igl::readOFF(argv[1], V, F);
	//	callback_load_mesh(viewer, argv[1]);
	}
	else
	{
		// Read mesh
		//callback_load_mesh(viewer, "../data/cube.off");
	}

	callback_key_down(viewer, '1', 0);

	//viewer.core.align_camera_center(V);
	viewer.launch();
}
