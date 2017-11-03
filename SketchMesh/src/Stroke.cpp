#include "Stroke.h"
#include <igl/unproject_onto_mesh.h>
#include <igl/unproject.h>
using namespace igl;
using namespace std;

std::chrono::steady_clock::time_point _time2;// = std::chrono::high_resolution_clock::now();
std::chrono::steady_clock::time_point _time1;// = std::chrono::high_resolution_clock::now();


Stroke::Stroke(const Eigen::MatrixXd &V_, const Eigen::MatrixXi &F_, igl::viewer::Viewer &v) :
	V(V_),
	F(F_),
	viewer(v) {
	stroke2DPoints = Eigen::MatrixX2d::Zero(1, 2);
	stroke3DPoints = Eigen::MatrixX3d::Zero(1, 3);
	_time1 = std::chrono::high_resolution_clock::now();

}

Stroke::~Stroke() {}

//To be used exclusively for drawing strokes for INITIAL MESHES (no extrusion)
void Stroke::strokeAddSegment(int mouse_x, int mouse_y) {
	//OpenGL has origin at left bottom, window(s) has origin at left top
	double x = mouse_x;
	double y = viewer.core.viewport(3) - mouse_y;
	if (!empty2D() && x == stroke2DPoints(stroke2DPoints.rows() - 1, 0) && y == stroke2DPoints(stroke2DPoints.rows() - 1, 1)) { //Check that the point is new compared to last time
		return;
	}

	if (!empty2D()) {
		_time2 = std::chrono::high_resolution_clock::now();
		auto timePast = std::chrono::duration_cast<std::chrono::nanoseconds>(_time2 - _time1).count();
		if (timePast < 100000) { //Don't add another segment before x nanoseconds
			return;
		}
	}

	Eigen::Matrix4f modelview = viewer.core.view * viewer.core.model;
	Eigen::RowVector3d pt(0, 0, 0);
	int faceID = -1;

	pt = igl::unproject(Eigen::Vector3f(x, y, 0.0f), modelview, viewer.core.proj, viewer.core.viewport).transpose().cast<double>();
	if (stroke2DPoints.rows() == 1 && empty2D()) {
		stroke2DPoints.row(0) << x, y;
		stroke3DPoints.row(0) << pt[0], pt[1], pt[2];
	} else {
		stroke2DPoints.conservativeResize(stroke2DPoints.rows() + 1, stroke2DPoints.cols());
		stroke2DPoints.row(stroke2DPoints.rows() - 1) << x, y;

		stroke3DPoints.conservativeResize(stroke3DPoints.rows() + 1, stroke3DPoints.cols());
		stroke3DPoints.row(stroke3DPoints.rows() - 1) << pt[0], pt[1], pt[2];
	}
	//using set_stroke_points will remove all previous strokes, using add_stroke_points might create duplicates
	viewer.data.set_stroke_points(stroke3DPoints);

	_time1 = std::chrono::high_resolution_clock::now(); //restart the "start" timer
}


//For drawing extrusion strokes. Need to start on the existing mesh
void Stroke::strokeAddSegmentExtrusion(int mouse_x, int mouse_y) {
	//OpenGL has origin at left bottom, window(s) has origin at left top
	double x = mouse_x;
	double y = viewer.core.viewport(3) - mouse_y;
	if (!empty2D() && x == stroke2DPoints(stroke2DPoints.rows() - 1, 0) && y == stroke2DPoints(stroke2DPoints.rows() - 1, 1)) { //Check that the point is new compared to last time
		return;
	}

	if (!empty2D()) {
		_time2 = std::chrono::high_resolution_clock::now();
		auto timePast = std::chrono::duration_cast<std::chrono::nanoseconds>(_time2 - _time1).count();
		if (timePast < 100000) { //Don't add another segment before x nanoseconds
			return;
		}
	}

	Eigen::Matrix4f modelview = viewer.core.view * viewer.core.model;
	Eigen::RowVector3d pt(0, 0, 0);
	int faceID = -1;

	if (dep < 0) { //first time
		Eigen::Vector3f bc;
		if (igl::unproject_onto_mesh(Eigen::Vector2f(x, y), modelview, viewer.core.proj, viewer.core.viewport, V, F, faceID, bc)) {
			pt = V.row(F(faceID, 0))*bc(0) + V.row(F(faceID, 1))*bc(1) + V.row(F(faceID, 2))*bc(2);
			Eigen::Vector3f proj = igl::project(pt.transpose().cast<float>().eval(), modelview, viewer.core.proj, viewer.core.viewport);
			dep = proj[2]; //0.85 is at center of meshes?
		}
	}
	pt = igl::unproject(Eigen::Vector3f(x, y, 0.95*dep), modelview, viewer.core.proj, viewer.core.viewport).transpose().cast<double>();

	if (stroke2DPoints.rows() == 1 && empty2D()) {
		stroke2DPoints.row(0) << x, y;
		stroke3DPoints.row(0) << pt[0], pt[1], pt[2];
	} else {
		stroke2DPoints.conservativeResize(stroke2DPoints.rows() + 1, stroke2DPoints.cols());
		stroke2DPoints.row(stroke2DPoints.rows() - 1) << x, y;

		stroke3DPoints.conservativeResize(stroke3DPoints.rows() + 1, stroke3DPoints.cols());
		stroke3DPoints.row(stroke3DPoints.rows() - 1) << pt[0], pt[1], pt[2];
	}
	//using set_stroke_points will remove all previous strokes, using add_stroke_points might create duplicates
	viewer.data.set_stroke_points(stroke3DPoints);

	_time1 = std::chrono::high_resolution_clock::now(); //restart the "start" timer
}

void Stroke::strokeReset() {
	stroke2DPoints.resize(1, 2);
	stroke2DPoints.setZero();
	stroke3DPoints.resize(1, 3);
	stroke3DPoints.setZero();
	dep = -1;
}

void Stroke::toLoop() {
	if(stroke2DPoints.rows() > 2) { //Don't do anything if we have only 1 line segment
		stroke2DPoints.conservativeResize(stroke2DPoints.rows() + 1, stroke2DPoints.cols());
		stroke2DPoints.row(stroke2DPoints.rows() - 1) << stroke2DPoints(0, 0), stroke2DPoints(0, 1);
		stroke3DPoints.conservativeResize(stroke3DPoints.rows() + 1, stroke3DPoints.cols());
		stroke3DPoints.row(stroke3DPoints.rows() - 1) << stroke3DPoints(0, 0), stroke3DPoints(0, 1), stroke3DPoints(0,2);
		//using set_stroke_points will remove all previous strokes, using add_stroke_points might create duplicates
		viewer.data.set_stroke_points(stroke3DPoints);
	}
}