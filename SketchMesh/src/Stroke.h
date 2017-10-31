#include <igl/viewer/Viewer.h>

class Stroke {
public:
	
	Stroke(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, const igl::viewer::Viewer &v);
	~Stroke();

private:
	const Eigen::MatrixXd &V;
	const Eigen::MatrixXi &F;
	const igl::viewer::Viewer &viewer;

	void strokeAddSegment(int mouse_x, int mouse_y);
	std::vector<Eigen::Matrix<double, 1, 3>> stroke3DPoints;
	std::vector<std::vector<int>> stroke2DPoints;
};