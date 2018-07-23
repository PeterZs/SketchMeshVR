#ifndef _Curve_Deformation_H_
#define _Curve_Deformation_H_
#include "Stroke.h"

class CurveDeformation {

public:
	static void startPullCurve(Stroke & _stroke, int _handle_ID, Eigen::MatrixXd & V, Eigen::MatrixXi & F);
	static void pullCurve(const Eigen::RowVector3d& pos, Eigen::MatrixXd& V);

	static bool smooth_deform_mode; //Used in the viewer menu to decide on smooth/sharp deformation
	static Eigen::MatrixXi EF, EV, FE;
	static std::vector<std::vector<int>> neighbors;

	struct PulledCurve {
		std::vector<int> vertices, fixed_vertices, edges, fixed_edges;
		Eigen::MatrixXi edge_triplets, vertex_triplets;


		 PulledCurve() {
			/*vertices.resize(0);
			fixed_vertices.resize(0);
			edges.resize(0);
			fixed_edges.resize(0);*/
			edge_triplets.resize(0, 3);
			vertex_triplets.resize(0, 3);
		}
	};
private:

	static double compute_curve_diag_length(Stroke & _stroke);
	static bool update_ROI(double drag_size);
	static void pullCurveTest(const Eigen::RowVector3d & pos, Eigen::MatrixXd & V, Eigen::VectorXi& edge_boundary_markers);
	static bool update_ROI_test(double drag_size, Eigen::MatrixXd & V, Eigen::VectorXi& edge_boundary_markers);
	static std::vector<int> collect_vertices_within_drag_length(double drag_size, Eigen::MatrixXd & V, Eigen::VectorXi& edge_boundary_markers);
	static void propagate(int prev_edge, int vert, double remaining_distance, std::vector<int>& vertices_in_range, Eigen::MatrixXd & V, Eigen::VectorXi& edge_boundary_markers);
	static void propagate(PulledCurve & curve, int prev_edge, int vert, int edge_marker, Eigen::VectorXi & edge_consumed, Eigen::VectorXi & edge_boundary_markers);
	static int find_next_edge(int vert, int prev_edge, int edge_marker, Eigen::VectorXi & edge_boundary_markers);
	static int find_edge(int start, int end);
	static void sort_by_distance(std::vector<int> vertices_in_range);
	static bool distance_sorter(int a, int b);
	static bool create_pulled_curve_by_propagation(int vert, Eigen::VectorXi & edge_consumed, PulledCurve & curve, Eigen::VectorXi& edge_boundary_markers);
	static std::vector<int> get_adjacent_seam_vertices(int vert, PulledCurve & curve);
	static int get_adjacent_seam_edge(int vert, int edge, std::vector<int> edges, Eigen::VectorXi & edge_boundary_markers);
	//static void compute_ROI_boundaries(int & ROI_1, int & ROI_2);
	//static void setup_fixed_indices(int ROI_1, int ROI_2);
	//static void setup_for_update_curve(Eigen::MatrixXd& V);
	static void setup_for_update_curve_test(PulledCurve & curve, Eigen::MatrixXd& V);
	static void setup_for_L1_position_step_test(PulledCurve & curve, Eigen::MatrixXd& V);
	//static void setup_for_L1_position_step(Eigen::MatrixXd& V);
	//static void update_curve(Eigen::MatrixXd & V);
	static void update_curve_test(PulledCurve & curve, Eigen::MatrixXd & V);
	static void solve_for_pos_and_rot_test(Eigen::MatrixXd & V, PulledCurve & curve);
	//static void solve_for_pos_and_rot(Eigen::MatrixXd& V);
	//static void update_rot();
	static void update_rot_test(PulledCurve & curve);
	static Eigen::Matrix3d compute_orthonormal(Eigen::Matrix3d & rot);
	//static void final_L1_pos(Eigen::MatrixXd & V);

	static void final_L1_pos_test(Eigen::MatrixXd & V, PulledCurve & curve);

	static Eigen::Matrix3d average_rot(Eigen::Matrix3d & r0, Eigen::Matrix3d & r1);

	
};
#endif

