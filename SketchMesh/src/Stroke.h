#ifndef _Stroke_
#define _Stroke_
#include <igl/viewer/Viewer.h>

class Stroke {
public:
	
	Stroke(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,  igl::viewer::Viewer &v);
	~Stroke();
	bool empty2D() const { return stroke2DPoints.isZero(); }
	bool empty3D() const { return stroke3DPoints.isZero(); }
	void strokeAddSegment(int mouse_x, int mouse_y);
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
	const Eigen::MatrixXd &V;
	const Eigen::MatrixXi &F;
	igl::viewer::Viewer &viewer;

	Eigen::MatrixX3d stroke3DPoints;
	Eigen::MatrixX2d stroke2DPoints;
	//std::vector<Eigen::Matrix<double, 1, 3>> stroke3DPoints;
	//std::vector<std::vector<int>> stroke2DPoints;
};

#endif