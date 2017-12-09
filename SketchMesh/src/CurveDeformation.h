#ifndef _Curve_Deformation_H_
#define _Curve_Deformation_H_
#include <Eigen/Core>
#include <Eigen/Geometry>
#include "Stroke.h"

class CurveDeformation {
public:
	static void startPullCurve(Stroke & _stroke, int handle_ID, int no_total_vertices);
	static void pullCurve(Eigen::RowVector3d& pos, Eigen::MatrixXd& V);

	static bool smooth_deform_mode;


private:
	static double current_max_drag_size;
	static double current_ROI_size;
	static double drag_size;
	static int handle_ID;
	static Eigen::RowVector3d start_pos;
	static Eigen::VectorXi fixed_indices;
	static Eigen::VectorXi fixed_indices_local;
	static int moving_vertex_ID;
	static double curve_diag_length;

	static double compute_curve_diag_length(Stroke & _stroke);
	static bool update_ROI(double drag_size);
	static void setup_for_update_curve(Eigen::MatrixXd& V);
	static void setup_for_L1_position_step(Eigen::MatrixXd& V);
	static void update_curve(Eigen::MatrixXd & V);
	static void solve_for_pos_and_rot(Eigen::MatrixXd& V);
	static void update_rot();
	static Eigen::Matrix3d compute_orthonormal(Eigen::Matrix3d & rot);
	static void final_L1_pos(Eigen::MatrixXd& V);
};
#endif

