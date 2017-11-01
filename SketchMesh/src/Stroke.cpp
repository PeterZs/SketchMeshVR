#include "Stroke.h"

using namespace igl;
using namespace std;

Stroke::Stroke(const Eigen::MatrixXd &V_, const Eigen::MatrixXi &F_,  igl::viewer::Viewer &v):
	V(V_),
	F(F_),
	viewer(v){
	stroke2DPoints = Eigen::MatrixX2d::Zero(1, 2);
	stroke3DPoints = Eigen::MatrixX3d::Zero(1, 3);
}

Stroke::~Stroke() {}

void Stroke::strokeAddSegment(int mouse_x, int mouse_y) {
	//OpenGL has origin at left bottom, window(s) has origin at left top
	double x = mouse_x;
	double y = viewer.core.viewport(3) - mouse_y;
	cout << x << " " << stroke2DPoints(stroke2DPoints.rows() - 1, 0) << " " << y << " " << stroke2DPoints(stroke2DPoints.rows() - 1, 0) << endl;
	cout << !empty2D() << endl;
	if (!empty2D() && x == stroke2DPoints(stroke2DPoints.rows()-1,0) && y == stroke2DPoints(stroke2DPoints.rows()-1,1)) { //Check that the point is new compared to last time
		return;
	}

	if (stroke2DPoints.rows()==1 && empty2D()) {
		cout << "empty" << endl;
		stroke2DPoints.row(0) << x, y;
		cout << stroke2DPoints.row(0) << endl;
	}
	else {

		//std::vector<int> newPoint;
		//newPoint.push_back(x);
		//newPoint.push_back(y);
		stroke2DPoints.conservativeResize(stroke2DPoints.rows()+1, stroke2DPoints.cols());
		//stroke2DPoints.row(stroke2DPoints.rows()-1) << x, y;
		cout << stroke2DPoints.rows() << " " << stroke2DPoints.cols() << endl;
		stroke2DPoints.row(stroke2DPoints.rows()-1) << x, y;
		//stroke2DPoints.push_back(newPoint);
		viewer.data.add_stroke_points(stroke2DPoints);
	}
}
