#ifndef _Curve_Deformation_H_
#define _Curve_Deformation_H_
#include <Eigen/Core>
#include <Eigen/Geometry>
#include "Stroke.h"

class CurveDeformation {
public:
	static void startPullCurve(Stroke& _stroke, int handle_ID);

	static void pullCurve(Eigen::RowVector3d& pos, Eigen::MatrixXd& V);


private:
	static double current_max_drag_size;
	static double CurveDeformation::current_ROI_size;
	static double CurveDeformation::drag_size;
	static Eigen::RowVector3d start_pos;
	static Eigen::VectorXi fixed_indices;
	static int moving_vertex_ID;
	static double curve_diag_length;

	static double compute_curve_diag_length(Stroke & _stroke);
	static bool update_ROI(double drag_size);
	static void setup_for_update_curve();
	static void setup_for_L1_position_step();
	static void update_curve();
};
#endif

