#include "Stroke.h"

using namespace igl;
using namespace std;

Stroke::Stroke(const Eigen::MatrixXd &V_, const Eigen::MatrixXi &F_, const igl::viewer::Viewer &v) :
	V(V_),
	F(F_),
	viewer(v) {}

Stroke::~Stroke() {}

void Stroke::strokeAddSegment(int mouse_x, int mouse_y) {
	//OpenGL has origin at left bottom, window(s) has origin at left top
	double x = mouse_x;
	double y = viewer.core.viewport(3) - mouse_y;

	std::vector<int> newPoint;
	newPoint.push_back(x);
	newPoint.push_back(y);
	stroke2DPoints.push_back(newPoint);
}
