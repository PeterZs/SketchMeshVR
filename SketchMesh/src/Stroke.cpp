#include "Stroke.h"
#include <igl/unproject_onto_mesh.h>
#include <igl/unproject.h>
#include <igl/triangle/triangulate.h>
#include <algorithm> 
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
	stroke_edges = Eigen::MatrixXi::Zero(0, 2);
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
	if (stroke2DPoints.rows() == 1 && empty2D()) { //Add first point
		stroke2DPoints.row(0) << x, y;
		stroke3DPoints.row(0) << pt[0], pt[1], pt[2];
	} else {
		stroke2DPoints.conservativeResize(stroke2DPoints.rows() + 1, stroke2DPoints.cols());
		stroke2DPoints.row(stroke2DPoints.rows() - 1) << x, y;

		stroke3DPoints.conservativeResize(stroke3DPoints.rows() + 1, stroke3DPoints.cols());
		stroke3DPoints.row(stroke3DPoints.rows() - 1) << pt[0], pt[1], pt[2];

		stroke_edges.conservativeResize(stroke_edges.rows() + 1, stroke_edges.cols());
		stroke_edges.row(stroke_edges.rows() - 1) << stroke2DPoints.rows() - 2, stroke2DPoints.rows() - 1;
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

	if (stroke2DPoints.rows() == 1 && empty2D()) { //Add first point
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

void Stroke::generateMeshFromStroke() {
	counter_clockwise(); //Ensure the stroke is counter-clockwise, handy later

	Eigen::MatrixX2d original_stroke2DPoints = stroke2DPoints;
	stroke2DPoints = resample_stroke(original_stroke2DPoints);


	Eigen::MatrixXd V2;
	Eigen::MatrixXi F2;
	igl::triangle::triangulate((Eigen::MatrixXd) stroke2DPoints, stroke_edges, Eigen::MatrixXd(0,0), "a0.005q", V2, F2);

}

void Stroke::counter_clockwise() {
	double total_area = 0;
	Eigen::Vector2d prev, next;
	prev = (Eigen::Vector2d) stroke2DPoints.row(stroke2DPoints.rows() - 1);
	for(int i = 0; i < stroke2DPoints.rows(); i++) {
		next = (Eigen::Vector2d) stroke2DPoints.row(i);
		total_area += (prev[1] + next[1]) * (next[0] - prev[0]);
		prev = next;
	}

	if(total_area > 0) { //reverse the vector
		stroke2DPoints.colwise().reverse();
	}

}

Eigen::Matrix2d Stroke::resample_stroke(Eigen::MatrixX2d & original_stroke2DPoints) {
	Eigen::MatrixX2d new_stroke2DPoints = Eigen::MatrixX2d::Zero(original_stroke2DPoints.rows(), 2);
	int nr_iterations = max(5, original_stroke2DPoints.rows() / 4);
	for(int i = 0; i < nr_iterations; i++) {
		move_to_middle(original_stroke2DPoints, new_stroke2DPoints);
		move_to_middle(new_stroke2DPoints, original_stroke2DPoints); //TODO: is this to save us a copy? basically performing an extra move_to_middle step
	}
	return new_stroke2DPoints;
}


void Stroke::move_to_middle(Eigen::MatrixX2d &positions, Eigen::MatrixX2d &new_positions) {
	int n = positions.rows();
	for(int i = 0; i < n; i++) {
		Eigen::Vector2d prev = positions.row((i - 1) % n);
		Eigen::Vector2d cur = positions.row(i%n);
		Eigen::Vector2d next = positions.row((i + 1) % n);

		new_positions(i, 0) = (cur[0] * 2 + prev[0] + next[0]) / 4;
		new_positions(i, 1) = (cur[1] * 2 + prev[1] + next[1]) / 4;
	}
}


double Stroke::total_stroke_length() {
	double total_length = 0;
	for(int i = 0; i < stroke2DPoints.rows(); i++) {
		total_length += (stroke2DPoints.row(i) - stroke2DPoints.row((i + 1) % stroke2DPoints.rows())).norm();
	}
	return total_length;
}