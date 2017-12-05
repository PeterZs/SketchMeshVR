#ifndef _Stroke_H_
#define _Stroke_H_
#include <igl/viewer/Viewer.h>
#include <vector>

class Stroke {
public:
	
	Stroke(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,  igl::viewer::Viewer &v);
	~Stroke();
	bool empty2D() const { return stroke2DPoints.isZero(); }
	bool empty3D() const { return stroke3DPoints.isZero(); }
	void strokeAddSegment(int mouse_x, int mouse_y);
	void strokeAddSegmentAdd(int mouse_x, int mouse_y);
	void strokeAddSegmentExtrusion(int mouse_x, int mouse_y);
	bool toLoop();
	void strokeReset();
	void generate3DMeshFromStroke(Eigen::VectorXi &vertex_boundary_markers);
	double total_stroke_length();
	Eigen::MatrixX3d get3DPoints();
	int get_vertex_idx_for_point(int i);
	int selectClosestVertex(int mouse_x, int mouse_y, double& closest_distance);
	double compute_stroke_diag();
	void update_Positions(Eigen::MatrixXd V);
	void snap_to_vertices();
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
	const Eigen::MatrixXd &V;
	const Eigen::MatrixXi &F;
	igl::viewer::Viewer &viewer;

	Eigen::MatrixX3d stroke3DPoints; //Used for screen output
	Eigen::MatrixX2d stroke2DPoints; //Used for early checking if point is new (in screen coordinates)
	Eigen::MatrixXi stroke_edges;
	double dep = -1;

	std::vector<int> closest_vert_bindings;

	void counter_clockwise();
	static Eigen::MatrixX2d resample_stroke(Eigen::MatrixX2d & original_stroke2DPoints);
	static void move_to_middle(Eigen::MatrixX2d &positions, Eigen::MatrixX2d &new_positions);
	static void generate_backfaces(Eigen::MatrixXi &faces, Eigen::MatrixXi &back_faces);
};

#endif