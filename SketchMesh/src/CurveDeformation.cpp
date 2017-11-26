#include <iostream>
#include <igl/unproject.h>

#include "CurveDeformation.h"

using namespace std;
using namespace igl;

double CurveDeformation::current_max_drag_size, CurveDeformation::current_ROI_size, CurveDeformation::drag_size, CurveDeformation::curve_diag_length;
Eigen::RowVector3d CurveDeformation::start_pos;
Eigen::VectorXi CurveDeformation::fixed_indices(0);
int CurveDeformation::moving_vertex_ID;
int no_vertices, no_ROI_vert = -1;
Eigen::SparseMatrix<double> A;
Eigen::SparseMatrix<double> A_L1;
Eigen::VectorXd B;
Eigen::SparseMatrix<double> B_L1;
Eigen::SparseLU<Eigen::SparseMatrix<double>> solverL1; //Solver for final vertex positions (with L1)


void CurveDeformation::startPullCurve(Stroke& _stroke, int handle_ID) {
	CurveDeformation::current_max_drag_size = -1.0;
	CurveDeformation::start_pos = (_stroke.get3DPoints()).row(handle_ID);
	CurveDeformation::moving_vertex_ID = handle_ID;
	CurveDeformation::current_ROI_size = 0.0;
	CurveDeformation::curve_diag_length = compute_curve_diag_length(_stroke);
	no_vertices = _stroke.get3DPoints().rows();
}

//pos is the unprojection from the position to where the user dragged the vertex
void CurveDeformation::pullCurve(Eigen::RowVector3d& pos, Eigen::MatrixXd& V) {
	Eigen::Vector3d move_vec = (pos - start_pos).transpose();

	drag_size = move_vec.norm();
	drag_size /= curve_diag_length;
	bool ROI_is_updated = false;
	if(!current_ROI_size || (drag_size - current_ROI_size > 0.01)) { //Take the deformation and roi size as percentages
		ROI_is_updated = update_ROI(drag_size);
	}

	if(no_ROI_vert == 1 || current_ROI_size == 0) { //If we have no other "free" vertices other than the handle vertex, we simply move it to the target position
		V.row(moving_vertex_ID) = pos;
	} else {
		if(ROI_is_updated) {
			setup_for_update_curve();
		}
		update_curve();
	}
}

double CurveDeformation::compute_curve_diag_length(Stroke& _stroke) {
	Eigen::Vector3d maxBB = _stroke.get3DPoints().colwise().maxCoeff();
	Eigen::Vector3d minBB = _stroke.get3DPoints().colwise().minCoeff();
	return (maxBB - minBB).norm();
}

bool CurveDeformation::update_ROI(double drag_size) {
	if(current_max_drag_size >= drag_size) {
		return false;
	}
	current_max_drag_size = drag_size;

	int no_ROI_vert_tmp = min(drag_size * no_vertices, ceil(((no_vertices - 1) / 2) - 1)); //Determine how many vertices to the left and to the right to have free (at most half-1 of all vertices on each side, always at least 1 vertex fixed)
	if(no_ROI_vert == no_ROI_vert_tmp) { //number of vertices in ROI didn't change
		return false;
	}

	current_ROI_size = drag_size;
	no_ROI_vert = no_ROI_vert_tmp;

	int ROI_1 = moving_vertex_ID - no_ROI_vert;
	int ROI_2 = (moving_vertex_ID + no_ROI_vert) % no_vertices;

	if(moving_vertex_ID < no_ROI_vert) { //ROI_1 will wrap around, manually perform modulo because negative modulo messes up
		ROI_1 = no_vertices - (no_ROI_vert - moving_vertex_ID);
	}

	vector<int> fixed;
	if(ROI_1 < ROI_2) {
		for(int i = 0; i <= ROI_1; i++) {
			fixed.push_back(i);
		}
		for(int i = ROI_2; i < no_vertices; i++) {
			fixed.push_back(i);
		}
	} else {
		for(int i = ROI_2; i <= ROI_1; i++) {
			fixed.push_back(i);
		}
	}

	fixed_indices = Eigen::VectorXi::Map(fixed.data(), fixed.size());

	return true;
}


void CurveDeformation::setup_for_update_curve() {
	A.resize(no_vertices + no_vertices * 3 + fixed_indices.size() + fixed_indices.size(), no_vertices + no_vertices);
	B.resize(no_vertices + no_vertices * 3 + fixed_indices.size() + fixed_indices.size());

	setup_for_L1_position_step();
}

void CurveDeformation::setup_for_L1_position_step() {
	int CONSTRAINT_WEIGHT = 10000;
	A_L1 = Eigen::SparseMatrix<double>(no_vertices + fixed_indices.size(), no_vertices);

	A_L1.insert(0, no_vertices - 1) = -0.5;
	A_L1.insert(0, 0) = 1;
	A_L1.insert(0, 1) = -0.5;
	for(int i = 1; i < no_vertices; i++) {
		A_L1.insert(i, i - 1) = -0.5;
		A_L1.insert(i, i) = 1;
		A_L1.insert(i, (i + 1)%no_vertices) = -0.5;
	}

	for(int i = 0; i < fixed_indices.size(); i++) {
		A_L1.insert(no_vertices + i, fixed_indices[i]) = CONSTRAINT_WEIGHT;
	}

	Eigen::SparseMatrix<double> A_L1_T = A_L1.transpose();
	solverL1.compute(A_L1_T*A_L1);
}

void CurveDeformation::update_curve() {

}