#ifndef _Stroke_
#define _Stroke_
#include <igl/viewer/Viewer.h>

class Stroke {
public:
	
	Stroke(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, const igl::viewer::Viewer &v);
	~Stroke();
	bool empty2D() const { return stroke2DPoints.empty(); }
	bool empty3D() const { return stroke3DPoints.empty(); }
	void strokeAddSegment(int mouse_x, int mouse_y);


private:
	const Eigen::MatrixXd &V;
	const Eigen::MatrixXi &F;
	const igl::viewer::Viewer &viewer;

	std::vector<Eigen::Matrix<double, 1, 3>> stroke3DPoints;
	std::vector<std::vector<int>> stroke2DPoints;
};

#endif