#include "Stroke.h"
#include "CleanStroke3D.h"
#include <igl/unproject_onto_mesh.h>
#include <igl/unproject.h>
#include <igl/project.h>
#include <igl/triangle/triangulate.h>
#include <igl/adjacency_list.h>
#include <igl/edge_topology.h>
#include <algorithm> 
#include <unordered_map>

#include <igl/per_vertex_normals.h>
#include <igl/slice.h>
#include "Plane.h"
#include <ctime>
using namespace igl;
using namespace std;

double min_inter_point_distance = 0.0001125;
int MAX_NR_TRIANGLES = 10000; //For the entire mesh (front + backside)
Eigen::RowVector3d red_color(1, 0, 0);
Eigen::RowVector3d blue_color(0, 0, 1);

Stroke::Stroke(Eigen::MatrixXd *V_, Eigen::MatrixXi* F_, int stroke_ID_) :
	V(V_),
	F(F_),
	stroke_ID(stroke_ID_) {
	stroke2DPoints = Eigen::MatrixXd::Zero(1, 2);
	stroke3DPoints = Eigen::MatrixX3d::Zero(1, 3);
	stroke3DPointsBack = Eigen::MatrixX3d::Zero(1, 3);
	stroke_edges = Eigen::MatrixXi(0, 2);
	faces_hit = Eigen::MatrixXi::Zero(1, 2);
	hand_pos_at_draw = Eigen::MatrixXd::Zero(1, 3);
	dep = Eigen::VectorXd::Zero(0);
	closest_vert_bindings.clear();
	has_points_on_mesh = false;
	has_been_outside_mesh = false;
	starts_on_mesh = true;
	ends_on_mesh = false;
	stroke_color = Eigen::RowVector3d(0.8*(rand() / (double)RAND_MAX), 0.8*(rand() / (double)RAND_MAX), 0.8*(rand() / (double)RAND_MAX));
	has_been_reversed = false;
}

Stroke::Stroke(const Stroke& origin) :
	V(origin.V),
	F(origin.F),
	stroke_ID(origin.stroke_ID),
	stroke2DPoints(origin.stroke2DPoints),
	stroke3DPoints(origin.stroke3DPoints),
	stroke3DPointsBack(origin.stroke3DPointsBack),
	stroke_edges(origin.stroke_edges),
	faces_hit(origin.faces_hit),
	hand_pos_at_draw(origin.hand_pos_at_draw),
	dep(origin.dep),
	closest_vert_bindings(origin.closest_vert_bindings),
	has_points_on_mesh(origin.has_points_on_mesh),
	starts_on_mesh(origin.starts_on_mesh),
	ends_on_mesh(origin.ends_on_mesh),
	has_been_outside_mesh(origin.has_been_outside_mesh),
	has_been_reversed(origin.has_been_reversed),
	stroke_color(origin.stroke_color),
	is_loop(origin.is_loop) {
}

Stroke& Stroke::operator=(Stroke other) {
	swap(other);
	return *this;
}

void Stroke::swap(Stroke & tmp) {//The pointers to V and F will always be the same for all stroke instances, so no need to copy
	std::swap(this->stroke_ID, tmp.stroke_ID);
	std::swap(this->stroke2DPoints, tmp.stroke2DPoints);
	std::swap(this->stroke3DPoints, tmp.stroke3DPoints);
	std::swap(this->stroke3DPointsBack, tmp.stroke3DPointsBack);
	std::swap(this->stroke_edges, tmp.stroke_edges);
	std::swap(this->faces_hit, tmp.faces_hit);
	std::swap(this->hand_pos_at_draw, tmp.hand_pos_at_draw);
	std::swap(this->closest_vert_bindings, tmp.closest_vert_bindings);
	std::swap(this->has_points_on_mesh, tmp.has_points_on_mesh);
	std::swap(this->starts_on_mesh, tmp.starts_on_mesh);
	std::swap(this->ends_on_mesh, tmp.ends_on_mesh);
	std::swap(this->has_been_outside_mesh, tmp.has_been_outside_mesh);
	std::swap(this->has_been_reversed, tmp.has_been_reversed);
	std::swap(this->stroke_color, tmp.stroke_color);
	std::swap(this->dep, tmp.dep);
	std::swap(this->is_loop, tmp.is_loop);
}

Stroke::Stroke() : V(&Eigen::MatrixXd(0, 3)), F(&Eigen::MatrixXi(0, 3)) {}
Stroke::~Stroke() {}

/** Used for DRAW. Will add a new 3D point to the stroke (if it is new compared to the last point, and didn't follow up too soon). **/
bool Stroke::addSegment(Eigen::Vector3f& pos, igl::opengl::glfw::Viewer &viewer) {
	if (!stroke3DPoints.isZero()) {
		if ((stroke3DPoints.row(stroke3DPoints.rows() - 1) - pos.transpose().cast<double>()).squaredNorm() < min_inter_point_distance) {
			return false;
		}

		if (stroke3DPoints.rows() > 10) {
			if ((stroke3DPoints.row(0) - pos.transpose().cast<double>()).squaredNorm() < (stroke3DPoints.row(0) - stroke3DPoints.row(1)).squaredNorm()*4.0) { //User is too close to beginning point. Don't sample, and block DRAW mode
				return true; //Block DRAW mode till the buttons are released again
			}
		}
	}

	if (stroke3DPoints.rows() == 1 && stroke3DPoints.isZero()) {
		stroke3DPoints.row(0) << pos[0], pos[1], pos[2];
	}
	else {
		stroke3DPoints.conservativeResize(stroke3DPoints.rows() + 1, Eigen::NoChange);
		stroke3DPoints.row(stroke3DPoints.rows() - 1) << pos[0], pos[1], pos[2];

		stroke_edges.conservativeResize(stroke_edges.rows() + 1, Eigen::NoChange);
		stroke_edges.bottomRows(1) << stroke3DPoints.rows() - 2, stroke3DPoints.rows() - 1;
	}
	closest_vert_bindings.push_back(stroke3DPoints.rows() - 1); //In the case of DRAW this will match the vertex indices, since we start from 0
	viewer.data().add_edges(stroke3DPoints.topRows(stroke3DPoints.rows() - 1), stroke3DPoints.middleRows(1, stroke3DPoints.rows() - 1), blue_color);
	if (stroke3DPoints.rows() > 10 && (stroke3DPoints.bottomRows(1) - stroke3DPoints.row(0)).squaredNorm() < (stroke3DPoints.row(1) - stroke3DPoints.row(0)).squaredNorm()*8.0) { //Show the closing line if the current point is close enough the first point (and we have already at least 10 samples)
		viewer.data().add_edges(stroke3DPoints.bottomRows(1), stroke3DPoints.row(0), blue_color);
	}
	return false;
}

/** Used for ADD. Will add a new 3D point to the stroke (if it is new compared to the last point, and didn't follow up too soon). **/
void Stroke::addSegmentAdd(Eigen::Vector3f& pos, igl::opengl::glfw::Viewer &viewer) {
	vector<igl::Hit> hits;

	if (igl::ray_mesh_intersect(pos, viewer.oculusVR.get_right_touch_direction(), *V, *F, hits)) { //Intersect the ray from the Touch controller with the mesh to get the 3D point
		if (hits.size() < 2) { //Hand was inside mesh during drawing
			return;
		}
		Eigen::Vector3d hit_pos = (*V).row((*F)(hits[0].id, 0))*(1.0 - hits[0].u - hits[0].v) + (*V).row((*F)(hits[0].id, 1))*hits[0].u + (*V).row((*F)(hits[0].id, 2))*hits[0].v;
		Eigen::Vector3d hit_pos_back = (*V).row((*F)(hits[1].id, 0))*(1.0 - hits[1].u - hits[1].v) + (*V).row((*F)(hits[1].id, 1))*hits[1].u + (*V).row((*F)(hits[1].id, 2))*hits[1].v;

		if (!stroke3DPoints.isZero() && hit_pos[0] == stroke3DPoints(stroke3DPoints.rows() - 1, 0) && hit_pos[1] == stroke3DPoints(stroke3DPoints.rows() - 1, 1) && hit_pos[2] == stroke3DPoints(stroke3DPoints.rows() - 1, 2)) {//Check that the point is new compared to last time
			return;
		}
		if ((stroke3DPoints.row(stroke3DPoints.rows() - 1) - hit_pos.transpose()).squaredNorm() < min_inter_point_distance) {
			return;
		}

		has_points_on_mesh = true;
		Eigen::Matrix4f modelview = viewer.oculusVR.get_start_action_view() * viewer.core.get_model();
		Eigen::Vector3f hit_pos_tmp = hit_pos.cast<float>();
		Eigen::Vector3d hit_pos2D = igl::project(hit_pos_tmp, modelview, viewer.core.get_proj(), viewer.core.viewport).cast<double>();

		if (stroke3DPoints.rows() == 1 && stroke3DPoints.isZero()) {
			stroke2DPoints.row(0) << hit_pos2D[0], hit_pos2D[1];
			stroke3DPoints.row(0) << hit_pos[0], hit_pos[1], hit_pos[2];
			stroke3DPointsBack.row(0) << hit_pos_back[0], hit_pos_back[1], hit_pos_back[2];
			faces_hit.row(0) << hits[0].id, hits[1].id;
			hand_pos_at_draw.row(0) << pos.transpose().cast<double>();
		}
		else {
			stroke2DPoints.conservativeResize(stroke2DPoints.rows() + 1, Eigen::NoChange);
			stroke2DPoints.row(stroke2DPoints.rows() - 1) << hit_pos2D[0], hit_pos2D[1];

			stroke3DPoints.conservativeResize(stroke3DPoints.rows() + 1, Eigen::NoChange);
			stroke3DPoints.row(stroke3DPoints.rows() - 1) << hit_pos[0], hit_pos[1], hit_pos[2];

			stroke3DPointsBack.conservativeResize(stroke3DPointsBack.rows() + 1, Eigen::NoChange);
			stroke3DPointsBack.row(stroke3DPointsBack.rows() - 1) << hit_pos_back[0], hit_pos_back[1], hit_pos_back[2];

			faces_hit.conservativeResize(faces_hit.rows() + 1, Eigen::NoChange);
			faces_hit.row(faces_hit.rows() - 1) << hits[0].id, hits[1].id;

			hand_pos_at_draw.conservativeResize(hand_pos_at_draw.rows() + 1, Eigen::NoChange);
			hand_pos_at_draw.row(hand_pos_at_draw.rows() - 1) << pos.transpose().cast<double>();
		}

		if (stroke3DPoints.rows() > 1) {
			viewer.data().add_edges(stroke3DPoints.block(stroke3DPoints.rows() - 2, 0, 1, 3), stroke3DPoints.block(stroke3DPoints.rows() - 1, 0, 1, 3), Eigen::RowVector3d(0, 1, 0));
		}
		ends_on_mesh = true;
		prev_point_was_on_mesh = true;
	}
	else {
		if (stroke3DPoints.rows() == 1) { //If we haven't any points yet, store the controller position & direction data so we can later create a looped cutting path
			pos_before_cut = pos.cast<double>();
			dir_before_cut = viewer.oculusVR.get_right_touch_direction().cast<double>();
			starts_on_mesh = false;
		}
		else if (prev_point_was_on_mesh) { //Point is outside off the mesh, after we have previously added on-mesh points. Store it for later creating a looped cutting path
			pos_after_cut = pos.cast<double>();
			dir_after_cut = viewer.oculusVR.get_right_touch_direction().cast<double>();
			prev_point_was_on_mesh = false;
		}
		else {
			has_been_outside_mesh = true; //We have an off-mesh point, that is not the start or end point
		}
		ends_on_mesh = false;
	}

	return;
}

/** Used for CUT. Will add a new 3D point to the stroke (if it is new compared to the last point, and didn't follow up too soon) and will also add its projection on the mesh as a 2D point. The indices of the front and backside faces that are hit will also be stored.
If the ray doesn't intersect the mesh, we will store its origin and direction and later possible use it to compute the start or end point. **/
void Stroke::addSegmentCut(Eigen::Vector3f& pos, igl::opengl::glfw::Viewer &viewer) {

	vector<igl::Hit> hits;

	if (igl::ray_mesh_intersect(pos, viewer.oculusVR.get_right_touch_direction(), (*V), (*F), hits)) { //Intersect the ray from the Touch controller with the mesh to get the 3D point
		if (hits.size() < 2) { //User had hand inside or behind mesh while cutting
			return;
		}

		Eigen::Vector3d hit_pos = (*V).row((*F)(hits[0].id, 0))*(1.0 - hits[0].u - hits[0].v) + (*V).row((*F)(hits[0].id, 1))*hits[0].u + (*V).row((*F)(hits[0].id, 2))*hits[0].v;
		Eigen::Vector3d hit_pos_back = (*V).row((*F)(hits[1].id, 0))*(1.0 - hits[1].u - hits[1].v) + (*V).row((*F)(hits[1].id, 1))*hits[1].u + (*V).row((*F)(hits[1].id, 2))*hits[1].v;

		//Early return (no point added)
		if (!stroke3DPoints.isZero() && hit_pos[0] == stroke3DPoints(stroke3DPoints.rows() - 1, 0) && hit_pos[1] == stroke3DPoints(stroke3DPoints.rows() - 1, 1) && hit_pos[2] == stroke3DPoints(stroke3DPoints.rows() - 1, 2)) {//Check that the point is new compared to last time
			return;
		}
		if (stroke2DPoints.rows() == 1 && dir_before_cut.isZero()) { //There is a intersecting stroke point before the vars for the 0th (off-mesh) have been set
			std::cerr << "Please draw the first point outside of the mesh" << endl;
			return;
		}
		if ((stroke3DPoints.row(stroke3DPoints.rows() - 1) - hit_pos.transpose()).squaredNorm() < min_inter_point_distance) {
			return;
		}

		has_points_on_mesh = true;
		Eigen::Matrix4f modelview = viewer.oculusVR.get_start_action_view() * viewer.core.get_model();
		Eigen::Vector3f hit_pos_tmp = hit_pos.cast<float>();
		Eigen::Vector3d hit_pos2D = igl::project(hit_pos_tmp, modelview, viewer.core.get_proj(), viewer.core.viewport).cast<double>();

		if (stroke2DPoints.rows() == 1 && !dir_before_cut.isZero() && stroke2DPoints.isZero()) {
			stroke2DPoints.row(0) << hit_pos2D[0], hit_pos2D[1];
			stroke3DPoints.row(0) << hit_pos[0], hit_pos[1], hit_pos[2];
			stroke3DPointsBack.row(0) << hit_pos_back[0], hit_pos_back[1], hit_pos_back[2];
			faces_hit.row(0) << hits[0].id, hits[1].id;
			stroke_color = Eigen::RowVector3d(1, 0, 0);
		}
		else {
			stroke2DPoints.conservativeResize(stroke2DPoints.rows() + 1, Eigen::NoChange);
			stroke2DPoints.row(stroke2DPoints.rows() - 1) << hit_pos2D[0], hit_pos2D[1];

			stroke3DPoints.conservativeResize(stroke3DPoints.rows() + 1, Eigen::NoChange);
			stroke3DPoints.row(stroke3DPoints.rows() - 1) << hit_pos[0], hit_pos[1], hit_pos[2];

			stroke3DPointsBack.conservativeResize(stroke3DPointsBack.rows() + 1, Eigen::NoChange);
			stroke3DPointsBack.row(stroke3DPointsBack.rows() - 1) << hit_pos_back[0], hit_pos_back[1], hit_pos_back[2];

			faces_hit.conservativeResize(faces_hit.rows() + 1, Eigen::NoChange);
			faces_hit.row(faces_hit.rows() - 1) << hits[0].id, hits[1].id;
		}

		if (stroke3DPoints.rows() > 2) {
			viewer.data().add_edges(stroke3DPoints.block(stroke3DPoints.rows() - 2, 0, 1, 3), stroke3DPoints.block(stroke3DPoints.rows() - 1, 0, 1, 3), Eigen::RowVector3d(1, 0, 1));
		}

		prev_point_was_on_mesh = true;
	}
	else if (stroke2DPoints.rows() > 1 && prev_point_was_on_mesh) { //We need to add the final point of the stroke, even though it is off the mesh (needed in order to wrap around to the backside)
		pos_after_cut = pos.cast<double>();
		dir_after_cut = viewer.oculusVR.get_right_touch_direction().cast<double>();
		prev_point_was_on_mesh = false;
	}
	else if (stroke2DPoints.rows() == 1) { //Add the first point, which should be outside the mesh. Refresh the "first" point until we get to the last one before we enter the mesh
		pos_before_cut = pos.cast<double>();
		dir_before_cut = viewer.oculusVR.get_right_touch_direction().cast<double>();
	}

	return;
}

/** Used for EXTRUDE. Extrusion base strokes need to be drawn entirely on the mesh (points outside of it will be ignored) and needs to surround at least one vertex.
Will add a new 3D point to the stroke (if it is new compared to the last point, and didn't follow up too soon) and will also add its projection on the mesh as a 2D point. Will also store the indices of the faces that are hit. The closest vertex bindings are handled in SurfacePath. **/
void Stroke::addSegmentExtrusionBase(Eigen::Vector3f& pos, igl::opengl::glfw::Viewer &viewer) {
	vector<igl::Hit> hits;

	if (igl::ray_mesh_intersect(pos, viewer.oculusVR.get_right_touch_direction(), (*V), (*F), hits)) { //Intersect the ray from the Touch controller with the mesh to get the 3D point
		if (hits.size() < 2) { //User had hand inside or behind mesh while drawing extrusion base stroke
			return;
		}
		Eigen::Vector3d hit_pos = (*V).row((*F)(hits[0].id, 0))*(1.0 - hits[0].u - hits[0].v) + (*V).row((*F)(hits[0].id, 1))*hits[0].u + (*V).row((*F)(hits[0].id, 2))*hits[0].v;


		if (!stroke3DPoints.isZero() && hit_pos[0] == stroke3DPoints(stroke3DPoints.rows() - 1, 0) && hit_pos[1] == stroke3DPoints(stroke3DPoints.rows() - 1, 1) && hit_pos[2] == stroke3DPoints(stroke3DPoints.rows() - 1, 2)) {//Check that the point is new compared to last time
			return;
		}
		else if ((stroke3DPoints.row(stroke3DPoints.rows() - 1) - hit_pos.transpose()).squaredNorm() < 0.25*min_inter_point_distance) {
			return;
		}

		has_points_on_mesh = true;
		Eigen::Matrix4f modelview = viewer.oculusVR.get_start_action_view() * viewer.core.get_model();
		Eigen::Vector3f hit_pos_tmp = hit_pos.cast<float>();
		Eigen::Vector3d hit_pos2D = igl::project(hit_pos_tmp, modelview, viewer.core.get_proj(), viewer.core.viewport).cast<double>();

		if (stroke2DPoints.rows() == 1 && stroke2DPoints.isZero()) {
			stroke2DPoints.row(0) << hit_pos2D[0], hit_pos2D[1];
			stroke3DPoints.row(0) << hit_pos[0], hit_pos[1], hit_pos[2];
			faces_hit.row(0) << hits[0].id, hits[1].id;
			hand_pos_at_draw.row(0) << pos.transpose().cast<double>();
			stroke_color = Eigen::RowVector3d(1, 0, 0);
		}
		else {
			stroke2DPoints.conservativeResize(stroke2DPoints.rows() + 1, Eigen::NoChange);
			stroke2DPoints.row(stroke2DPoints.rows() - 1) << hit_pos2D[0], hit_pos2D[1];

			stroke3DPoints.conservativeResize(stroke3DPoints.rows() + 1, Eigen::NoChange);
			stroke3DPoints.row(stroke3DPoints.rows() - 1) << hit_pos[0], hit_pos[1], hit_pos[2];

			faces_hit.conservativeResize(faces_hit.rows() + 1, Eigen::NoChange);
			faces_hit.row(faces_hit.rows() - 1) << hits[0].id, hits[1].id;

			hand_pos_at_draw.conservativeResize(hand_pos_at_draw.rows() + 1, Eigen::NoChange);
			hand_pos_at_draw.row(hand_pos_at_draw.rows() - 1) << pos.transpose().cast<double>();
		}
	}
	else {
		if (has_points_on_mesh) {
			has_been_outside_mesh = true;
		}
		return;
	}

	if (stroke3DPoints.rows() > 1) {
		viewer.data().add_edges(stroke3DPoints.block(stroke3DPoints.rows() - 2, 0, 1, 3), stroke3DPoints.block(stroke3DPoints.rows() - 1, 0, 1, 3), Eigen::RowVector3d(0, 0, 1));
	}
}

/** Used for EXTRUDE. Will add a new 3D point to the stroke (if it is new compared to the last point, and didn't follow up too soon) and will also add its projection as a 2D point. The closest vertex bindings are handled in SurfacePath. **/
void Stroke::addSegmentExtrusionSilhouette(Eigen::Vector3f& pos, igl::opengl::glfw::Viewer &viewer) {
	Eigen::RowVector3d pt2D;
	Eigen::Matrix4f modelview = viewer.oculusVR.get_start_action_view() * viewer.core.get_model();

	Eigen::RowVector3d pos_in = pos.cast<double>().transpose();
	igl::project(pos_in, modelview, viewer.core.get_proj(), viewer.core.viewport, pt2D);


	if (!stroke2DPoints.isZero() && pt2D[0] == stroke2DPoints(stroke2DPoints.rows() - 1, 0) && pt2D[1] == stroke2DPoints(stroke2DPoints.rows() - 1, 1)) {//Check that the point is new compared to last time
		return;
	}
	else if ((stroke3DPoints.row(stroke3DPoints.rows() - 1) - pos.transpose().cast<double>()).squaredNorm() < 0.25*min_inter_point_distance) {
		return;
	}

	if (stroke2DPoints.rows() == 1 && empty2D()) {
		stroke2DPoints.row(0) << pt2D[0], pt2D[1];
		stroke3DPoints.row(0) << pos[0], pos[1], pos[2];
	}
	else {
		stroke2DPoints.conservativeResize(stroke2DPoints.rows() + 1, Eigen::NoChange);
		stroke2DPoints.row(stroke2DPoints.rows() - 1) << pt2D[0], pt2D[1];

		stroke3DPoints.conservativeResize(stroke3DPoints.rows() + 1, Eigen::NoChange);
		stroke3DPoints.row(stroke3DPoints.rows() - 1) << pos[0], pos[1], pos[2];
	}

	if (stroke3DPoints.rows() > 1) {
		viewer.data().add_edges(stroke3DPoints.block(stroke3DPoints.rows() - 2, 0, 1, 3), stroke3DPoints.block(stroke3DPoints.rows() - 1, 0, 1, 3), Eigen::RowVector3d(0.2, 0.6, 1));
	}

}

/** Used for CUT. Adds the first point, which is the last point outside of the mesh before cut start. Doesn't get drawn. **/
void Stroke::prepend_first_point(igl::opengl::glfw::Viewer &viewer) {
	Eigen::Vector3d first_point = pos_before_cut + (stroke3DPoints.row(0).transpose() - pos_before_cut).dot(dir_before_cut.normalized()) * dir_before_cut.normalized(); //The closest point Pr along a line that starts from P1 and does in direction dir to point P2 is as follows: Pr = P1 + (P2 - P1).dot(dir) * dir with dir normalized
	Eigen::Matrix4f modelview = viewer.oculusVR.get_start_action_view() * viewer.core.get_model();
	Eigen::Vector3f first_point_tmp = first_point.cast<float>();
	Eigen::Vector3d hit_pos2D = igl::project(first_point_tmp, modelview, viewer.core.get_proj(), viewer.core.viewport).cast<double>();

	stroke2DPoints.conservativeResize(stroke2DPoints.rows() + 1, Eigen::NoChange);
	Eigen::MatrixXd old = stroke2DPoints.block(0, 0, stroke2DPoints.rows() - 1, 2);
	stroke2DPoints.block(1, 0, stroke2DPoints.rows() - 1, 2) = old;
	stroke2DPoints.row(0) << hit_pos2D[0], hit_pos2D[1];

	stroke3DPoints.conservativeResize(stroke3DPoints.rows() + 1, Eigen::NoChange);
	old = stroke3DPoints.block(0, 0, stroke3DPoints.rows() - 1, 3);
	stroke3DPoints.block(1, 0, stroke3DPoints.rows() - 1, 3) = old;
	stroke3DPoints.row(0) = first_point;

	faces_hit.conservativeResize(faces_hit.rows() + 1, Eigen::NoChange);
	Eigen::MatrixXi old_faces_hit = faces_hit.block(0, 0, faces_hit.rows() - 1, 2);
	faces_hit.block(1, 0, faces_hit.rows() - 1, 2) = old_faces_hit;
	faces_hit.row(0) << -1, -1;
}

/** Used for CUT. Adds the final point, which is the first point outside of the mesh. Doesn't get drawn. **/
void Stroke::append_final_point(igl::opengl::glfw::Viewer &viewer) {
	Eigen::Vector3d last_point = pos_after_cut + (stroke3DPoints.row(stroke3DPoints.rows() - 1).transpose() - pos_after_cut).dot(dir_after_cut.normalized()) * dir_after_cut.normalized(); //The closest point Pr along a line that starts from P1 and goes in direction dir to point P2 is as follows: Pr = P1 + (P2 - P1).dot(dir) * dir with dir normalized
	Eigen::Matrix4f modelview = viewer.oculusVR.get_start_action_view() * viewer.core.get_model();
	Eigen::Vector3f last_point_tmp = last_point.cast<float>();
	Eigen::Vector3d hit_pos2D = igl::project(last_point_tmp, modelview, viewer.core.get_proj(), viewer.core.viewport).cast<double>();

	stroke2DPoints.conservativeResize(stroke2DPoints.rows() + 1, Eigen::NoChange);
	stroke2DPoints.row(stroke2DPoints.rows() - 1) << hit_pos2D[0], hit_pos2D[1];

	stroke3DPoints.conservativeResize(stroke3DPoints.rows() + 1, Eigen::NoChange);
	stroke3DPoints.row(stroke3DPoints.rows() - 1) = last_point;

	faces_hit.conservativeResize(faces_hit.rows() + 1, Eigen::NoChange);
	faces_hit.row(faces_hit.rows() - 1) << -1, -1;
}

bool Stroke::toLoop() {
	if (stroke3DPoints.rows() > 2) { //Don't do anything if we have only 1 line segment
		stroke3DPoints.conservativeResize(stroke3DPoints.rows() + 1, Eigen::NoChange);
		stroke3DPoints.row(stroke3DPoints.rows() - 1) << stroke3DPoints.row(0);
		if (closest_vert_bindings.size() > 0) {
			closest_vert_bindings.push_back(closest_vert_bindings[0]);
		}
		if (stroke_edges.rows() > 0) {
			stroke_edges.conservativeResize(stroke_edges.rows() + 1, Eigen::NoChange);
			stroke_edges.bottomRows(1) << closest_vert_bindings[closest_vert_bindings.size() - 2], closest_vert_bindings.back();
		}
		if (hand_pos_at_draw.rows() > 1) {
			hand_pos_at_draw.conservativeResize(hand_pos_at_draw.rows() + 1, Eigen::NoChange);
			hand_pos_at_draw.bottomRows(1) << hand_pos_at_draw.row(0);
		}
		if (faces_hit.rows() > 1) {
			faces_hit.conservativeResize(faces_hit.rows() + 1, Eigen::NoChange);
			faces_hit.bottomRows(1) << faces_hit.row(0);
		}
		is_loop = true;
		return true;
	}
	is_loop = false;
	return false;
}

void Stroke::generate3DMeshFromStroke(Eigen::VectorXi &edge_boundary_markers, Eigen::VectorXi& vertex_is_fixed, Eigen::MatrixXd& mesh_V, Eigen::MatrixXi& mesh_F, igl::opengl::glfw::Viewer &viewer) {
	double mean_sample_dist = 0.0;
	for (int i = 0; i < stroke3DPoints.rows() - 1; i++) {
		mean_sample_dist += (stroke3DPoints.row(i) - stroke3DPoints.row((i + 1) % stroke3DPoints.rows())).norm();
	}
	mean_sample_dist /= stroke3DPoints.rows();
	Eigen::MatrixXd closing_stroke3DPoints(2, 3);
	closing_stroke3DPoints.row(0) = stroke3DPoints.row(stroke3DPoints.rows() - 2);
	closing_stroke3DPoints.row(1) = stroke3DPoints.row(0);
	Eigen::MatrixXd resampled_3DPoints = CleanStroke3D::resample_by_length_sub(closing_stroke3DPoints, 0, 1, mean_sample_dist);

	Eigen::MatrixXd new_3DPoints = stroke3DPoints.topRows(stroke3DPoints.rows() - 1);
	new_3DPoints.conservativeResize(new_3DPoints.rows() + resampled_3DPoints.rows(), Eigen::NoChange);
	new_3DPoints.bottomRows(resampled_3DPoints.rows()) = resampled_3DPoints;

	Eigen::MatrixXd V2_tmp, V2;
	Eigen::MatrixXi F2, F2_back, vertex_markers;
	Eigen::RowVector3d center;
	Eigen::Vector3d xvec, yvec;
	int count = 0;
	double resampling_distance = sqrt(min_inter_point_distance);
	do{
		count++;
		if (count > 1) { //Previous triangulation was too finegrained, have to increase resampling distance
			resampling_distance *= pow(1.15, 2*double(abs((MAX_NR_TRIANGLES/2.0)-F2.rows()) / (MAX_NR_TRIANGLES / 2.0))); //Multiply the exponent by a factor 2 to ensure faster convergence
		}
		std::cout << "Step0 " << resampling_distance << std::endl;
		Eigen::MatrixXd resampled_new_3DPoints = CleanStroke3D::resample_by_length_sub(new_3DPoints, 0, new_3DPoints.rows() - 1, resampling_distance);
		closest_vert_bindings.clear();

		for (int i = 0; i < resampled_new_3DPoints.rows(); i++) {
			closest_vert_bindings.push_back(i);
		}
		set3DPoints(resampled_new_3DPoints);

		center = stroke3DPoints.colwise().mean();
		stroke3DPoints = stroke3DPoints.rowwise() - center; //Zero mean
		Eigen::MatrixXd cov = (1.0 / (stroke3DPoints.rows() - 1)) * (stroke3DPoints.transpose()*stroke3DPoints);
		EigenSolver<MatrixXd> es;
		es.compute(cov);
		Eigen::VectorXd eigenvals = es.eigenvalues().real();

		std::vector<int> idx(eigenvals.rows());
		size_t n(0);
		std::generate(std::begin(idx), std::end(idx), [&] {return n++; }); //Fill the vector idx with increasing integers starting from 0
		std::sort(std::begin(idx), std::end(idx), [&](int i1, int i2) {return eigenvals[i1] > eigenvals[i2]; }); //Sort the vector of indices based on the descending values of eigenvalues

		xvec = es.eigenvectors().real().col(idx[0]);
		yvec = es.eigenvectors().real().col(idx[1]);
		Eigen::Vector3d zvec = es.eigenvectors().real().col(idx[2]);
		Eigen::MatrixXd new_axes(3, 3);
		new_axes << xvec, yvec, zvec;

		Eigen::MatrixX3d projected_points = stroke3DPoints * new_axes;
		dep = projected_points.col(2);

		stroke2DPoints = projected_points.leftCols(2);
		TaubinFairing2D(stroke2DPoints, 5);

		counter_clockwise();  //Ensure the stroke is counter-clockwise, handy later

		Eigen::MatrixXi edge_markers, stroke_edges(stroke2DPoints.rows(), 2);
		double mean_squared_sample_dist = 0.0;

		for (int i = 0; i < stroke2DPoints.rows(); i++) {
			stroke_edges.row(i) << i, ((i + 1) % stroke2DPoints.rows());
			mean_squared_sample_dist += (stroke2DPoints.row(i) - stroke2DPoints.row((i + 1) % stroke2DPoints.rows())).squaredNorm();
		}
		mean_squared_sample_dist /= stroke2DPoints.rows();

		std::cout << "Before triangulate" << std::endl;
		//Set a inter-sample distance dependent maximum triangle area. Setting this too small will result in inflation in the wrong direction.
		igl::triangle::triangulate((Eigen::MatrixXd) stroke2DPoints, stroke_edges, Eigen::MatrixXd(0, 0), Eigen::MatrixXi::Constant(stroke2DPoints.rows(), 1, 1), Eigen::MatrixXi::Constant(stroke_edges.rows(), 1, 1), "QYq25a" + to_string(0.5*(mean_squared_sample_dist)), V2_tmp, F2, vertex_markers, edge_markers);
		std::cout << " Done triangulating " << F2.rows() << " triangles." << std::endl;

	} while (F2.rows() > (MAX_NR_TRIANGLES/2));

	V2 = center.replicate(V2_tmp.rows(), 1);

	for (int i = 0; i < V2_tmp.rows(); i++) {
		V2.row(i) += xvec.transpose()*V2_tmp(i, 0);
		V2.row(i) += yvec.transpose()*V2_tmp(i, 1);
	}

	generate_backfaces(F2, F2_back);
	int nr_boundary_vertices = (vertex_markers.array() == 1).count();
	int original_size = V2.rows();
	V2.conservativeResize(2 * original_size - nr_boundary_vertices, Eigen::NoChange); //Increase size, such that boundary vertices only get included one
	unordered_map<int, int> backside_vertex_map;
	int count2 = 0;
	for (int i = 0; i < original_size; i++) {
		if (vertex_markers(i) != 1) { //If the vertex is NOT on the boundary, duplicate it
			V2.row(original_size + count2) << V2.row(i);
			backside_vertex_map.insert({ i, original_size + count2 });
			backside_vertex_map.insert({ original_size + count2, i });
			count2++;
		}
		else {
			backside_vertex_map.insert({ i, i });
		}
	}

	//Add backside faces, using original vertices for boundary and copied vertices for interior
	original_size = F2_back.rows();
	F2.conservativeResize(original_size * 2, Eigen::NoChange);
	for (int i = 0; i < original_size; i++) {
		for (int j = 0; j < 3; j++) {
			auto got = backside_vertex_map.find(F2_back(i, j));
			F2(original_size + i, j) = got->second;
		}
	}

	Eigen::MatrixXd N_Faces, N_Vertices;
	igl::per_vertex_normals(V2, F2, PER_VERTEX_NORMALS_WEIGHTING_TYPE_UNIFORM, N_Vertices);
	vertex_is_fixed.resize(V2.rows());
	for (int i = 0; i < V2.rows(); i++) {
		if (i >= vertex_markers.rows()) { //vertex can't be boundary (it's on backside)
			V2.row(i) = V2.row(i) + 0.05*N_Vertices.row(i);
			vertex_is_fixed[i] = 0;
		}
		else {
			if (vertex_markers(i) == 1) { //Don't change boundary vertices
				vertex_is_fixed[i] = 1;
				//Don't need to change stroke3DPoints here because they're the same as V2's points.
				continue;
			}
			V2.row(i) = V2.row(i) + 0.05*N_Vertices.row(i);
			vertex_is_fixed[i] = 0;
		}
	}

	mesh_V = V2;
	mesh_F = F2;

	if (!is_edge_manifold(*F)) {
		stroke3DPoints.conservativeResize(stroke3DPoints.rows() + 1, Eigen::NoChange);
		stroke3DPoints.bottomRows(1) = stroke3DPoints.row(0);
		stroke3DPoints.rowwise() += center;
		closest_vert_bindings.push_back(closest_vert_bindings[0]);

		return;
	}
	Eigen::MatrixXi EV, FE, EF;
	igl::edge_topology(*V, *F, EV, FE, EF);	

	edge_boundary_markers.conservativeResize(EV.rows());
	edge_boundary_markers.setZero();
	for (int i = 0; i < EV.rows(); i++) {
		if (vertex_is_fixed[EV(i, 0)] && vertex_is_fixed[EV(i, 1)]) { //If both ends of the edge are fixed, then the edge must be on the boundary
			edge_boundary_markers[i] = 1;
		}
	}

	//Make vertex bindings and 3DPoints looped again
	stroke3DPoints.conservativeResize(stroke3DPoints.rows() + 1, Eigen::NoChange);
	stroke3DPoints.bottomRows(1) = stroke3DPoints.row(0);
	stroke3DPoints.rowwise() += center;
	closest_vert_bindings.push_back(closest_vert_bindings[0]);

	return;
}

void Stroke::project_with_PCA_given_target(Eigen::Vector3d target_vec) {
	Eigen::MatrixX3d projected_points;
	Eigen::RowVector3d center = stroke3DPoints.colwise().mean();
	stroke3DPoints = stroke3DPoints.rowwise() - center; //Zero mean
	Eigen::MatrixXd cov = (1.0 / (stroke3DPoints.rows() - 1)) * (stroke3DPoints.transpose()*stroke3DPoints);
	EigenSolver<MatrixXd> es;
	es.compute(cov);
	Eigen::VectorXd eigenvals = es.eigenvalues().real();

	std::vector<int> idx(eigenvals.rows());
	size_t n(0);
	std::generate(std::begin(idx), std::end(idx), [&] {return n++; }); //Fill the vector idx with increasing integers starting from 0
	std::sort(std::begin(idx), std::end(idx), [&](int i1, int i2) {return eigenvals[i1] > eigenvals[i2]; }); //Sort the vector of indices based on the descending values of eigenvalues
	

	double proj, max_proj = 0.0;
	int max_idx;
	for (int i = 0; i < 3; i++) {
		proj = abs(target_vec.dot(es.eigenvectors().real().col(idx[i])));
		if (proj > max_proj) {
			max_proj = proj;
			max_idx = i;
		}
	}

	Eigen::Vector3d xvec = es.eigenvectors().real().col(idx[max_idx]);
	Eigen::Vector3d yvec = (max_idx == 0) ? es.eigenvectors().real().col(idx[1]) : es.eigenvectors().real().col(idx[0]);
	Eigen::Vector3d zvec = (max_idx == 0) ? es.eigenvectors().real().col(idx[2]) : es.eigenvectors().real().col(!(max_idx - 1) + 1);

	Eigen::MatrixXd new_axes(3, 3);
	new_axes << xvec, yvec, zvec;

	projected_points = stroke3DPoints * new_axes;

	for (int i = 0; i < stroke3DPoints.rows(); i++) {
		stroke3DPoints.row(i) = center;
		stroke3DPoints.row(i) += xvec.transpose()*projected_points(i, 0);
		stroke3DPoints.row(i) += yvec.transpose()*projected_points(i, 1);
	}
}

/** Makes the stroke counter clockwise. **/
void Stroke::counter_clockwise() {
	double total_area = 0;
	Eigen::Vector2d prev, next;
	prev = (Eigen::Vector2d) stroke2DPoints.row(stroke2DPoints.rows() - 1);
	for (int i = 0; i < stroke2DPoints.rows(); i++) {
		next = (Eigen::Vector2d) stroke2DPoints.row(i);
		total_area += (prev[1] + next[1]) * (next[0] - prev[0]);
		prev = next;
	}

	if (total_area > 0) { //reverse the vector
		stroke2DPoints = stroke2DPoints.colwise().reverse().eval();
		stroke3DPoints = stroke3DPoints.colwise().reverse().eval();
		if (hand_pos_at_draw.rows() > 1) {
			hand_pos_at_draw = hand_pos_at_draw.colwise().reverse().eval();
		}
		if (faces_hit.rows() > 1) {
			faces_hit = faces_hit.colwise().reverse().eval();
		}
		has_been_reversed = true;
	}
}

/** Taubin fairing will make sure that the shape doesn't shrink when smoothing the curve **/
void Stroke::TaubinFairing2D(Eigen::MatrixXd& original_stroke2DPoints, int n) {
	for (int i = 0; i < n; i++) {
		smooth_sub(original_stroke2DPoints, 0.63139836);
		smooth_sub(original_stroke2DPoints, -0.6739516);
	}
}

void Stroke::smooth_sub(Eigen::MatrixXd& points, double direction) {
	int n = points.rows();
	Eigen::MatrixXd new_positions(n, 2);
	Eigen::RowVector2d prev, cur, next;
	for (int i = 0; i < n; i++) {
		prev = points.row(((i - 1) + n) % n);
		cur = points.row(i);
		next = points.row(((i + 1) + n) % n);
		new_positions.row(i) = to_sum_of_vectors(cur, prev, next, direction);
	}

	for (int i = 0; i < n; i++) {
		points.row(i) = new_positions.row(i);
	}
}

Eigen::RowVector2d Stroke::to_sum_of_vectors(Eigen::RowVector2d vert, Eigen::RowVector2d prev, Eigen::RowVector2d next, double direction) {
	Eigen::RowVector2d total_vec(0, 0);
	double total_w = 0.0;
	double w = 0.0;
	Eigen::RowVector2d vec;

	vec = prev - vert;
	w = 1.0 / vec.norm();
	total_vec += vec*w;
	total_w += w;

	vec = next - vert;
	w = 1.0 / vec.norm();
	total_vec += vec*w;
	total_w += w;

	total_vec *= (1.0 / total_w);
	return (vert + total_vec*direction);
}

void Stroke::generate_backfaces(Eigen::MatrixXi &faces, Eigen::MatrixXi &back_faces) {
	back_faces = faces.rowwise().reverse().eval();
}

/** Returns the ID of the stroke's 3D point that is closest to where the user clicked. Stores the distance_to_vert to this point in closest_distance. If the user clicked too far away from any of the stroke's points, it will return -1. **/
int Stroke::selectClosestVertex(Eigen::Vector3f pos, double& closest_distance) {
	double closest_dist = INFINITY, dist;
	int closest_ID;

	for (int i = 0; i < closest_vert_bindings.size(); i++) {
		dist = ((*V).row(closest_vert_bindings[i]).transpose().cast<float>() - pos).squaredNorm();
		if (dist < closest_dist) {
			closest_dist = dist;
			closest_ID = i;
		}
	}

	closest_dist = sqrt(closest_dist);
	double stroke_diag = compute_stroke_diag();
	closest_distance = closest_dist;
	if (closest_vert_bindings.size() > 2 && closest_dist > 0.08 * stroke_diag) { //Don't do pulling when the user clicks too far away from any curve (treat as navigation movement)
		return -1;
	}
	else if (closest_vert_bindings.size() == 2 && closest_dist > 0.08 * ((*V).row(closest_vert_bindings[0]).norm() * 2)) { //The stroke consist of a single point (that's "looped") so the diag will be 0. Use the point's norm instead
		return -1;
	}
	return closest_vert_bindings[closest_ID];
}

double Stroke::compute_stroke_diag() {
	Eigen::MatrixXd stroke_points;
	Eigen::VectorXi row_slice = Map<Eigen::VectorXi >(closest_vert_bindings.data(), closest_vert_bindings.size());
	Eigen::VectorXi col_slice(3);
	col_slice.col(0) << 0, 1, 2;
	igl::slice((*V), row_slice, col_slice, stroke_points);

	Eigen::Vector3d maxBB = stroke_points.colwise().maxCoeff();
	Eigen::Vector3d minBB = stroke_points.colwise().minCoeff();
	return (maxBB - minBB).norm();
}

void Stroke::update_Positions(Eigen::MatrixXd V, bool structure_changed) {
	stroke3DPoints.resize(closest_vert_bindings.size(), Eigen::NoChange);
	for (int i = 0; i < closest_vert_bindings.size(); i++) { //Iterate over the (updated) closest_vert_bindings instead of over stroke3DPoints
		stroke3DPoints.row(i) = V.row(closest_vert_bindings[i]);
	}
	//In the case of extrusion silhouette strokes, closest_vert_bindings isn't looped. Don't make stroke3DPoints looped, because we already account for it not being a loop when drawing the curves

	if (structure_changed) {
		Eigen::MatrixXi EV, FE, EF;
		igl::edge_topology(V, *F, EV, FE, EF);

		int start, end, equal_pos;
		Eigen::VectorXi col1Equals, col2Equals;
		stroke_edges.resize(0, 2);
		for (int i = 0; i < closest_vert_bindings.size() - 1; i++) { //Closest_vert_bindings is already looped, so don't double loop it
			start = closest_vert_bindings[i];
			end = closest_vert_bindings[(i + 1) % closest_vert_bindings.size()];
			col1Equals = EV.col(0).cwiseEqual(min(start, end)).cast<int>();
			col2Equals = EV.col(1).cwiseEqual(max(start, end)).cast<int>();
			int maxval = (col1Equals + col2Equals).maxCoeff(&equal_pos); //Find the row that contains both vertices of this edge
			if (maxval == 2) {
				stroke_edges.conservativeResize(stroke_edges.rows() + 1, Eigen::NoChange);
				stroke_edges.bottomRows(1) << start, end;
			}
		}
		if (!is_loop) {
			stroke_edges.conservativeResize(stroke_edges.rows() + 1, Eigen::NoChange);
			stroke_edges.bottomRows(1) << closest_vert_bindings[closest_vert_bindings.size() - 2], closest_vert_bindings.back(); //Manually make it looped again, with a connection between the last 2 vertex bindings (even if there is not an actual edge connecting them, e.g. with interrupted strokes due to a cut). Drawing the strokes assumes they're always looped, and then takes the correct number of rows based on is_loop
		}
	}
}

/** Remaps the stroke's vertex bindings after the mesh topology has changed.
	Note that this method does not update the stroke3DPoints, use update_Positions() for that. **/
bool Stroke::update_vert_bindings(Eigen::VectorXi & new_mapped_indices, Eigen::MatrixXi& replacing_vertex_bindings) {
	vector<int> new_bindings;
	int first_included_after_remove = -1;
	bool points_were_removed = false, no_tracked_point_yet = true, stays_continuous = false, originally_is_loop = is_loop;

	for (int i = 0; i < replacing_vertex_bindings.rows(); i++) {
		if (replacing_vertex_bindings(i, 0) == stroke_ID) {
			auto loc = std::find(closest_vert_bindings.begin(), closest_vert_bindings.end(), replacing_vertex_bindings(i, 1));
			int idx = distance(closest_vert_bindings.begin(), loc);
			loc = std::find(closest_vert_bindings.begin(), closest_vert_bindings.end(), replacing_vertex_bindings(i, 2));
			int idx2 = distance(closest_vert_bindings.begin(), loc);
			if (abs(idx2 - idx) != 1) { //Take absolute in case vertices are found in reversed order (e.g. first 2 vertex bindings of the original stroke which has been cut, so the first stroke vertex and the intersection vertex to the cut surface). Find will return the position of the first found one 
				std::cerr << std::endl << replacing_vertex_bindings(i, 1) << "   " << replacing_vertex_bindings(i, 2) << " Something went wrong. Multiple vertices inserted on 1 edge?" << std::endl;
			}
			else {
				closest_vert_bindings.insert(loc, replacing_vertex_bindings(i, 3)); //Insert the new middle vertex between the 2 existing edge vertices
			}
		}
	}

	for (int i = 0; i < closest_vert_bindings.size() - 1; i++) {	//Closest_vert_bindings is always a loop
		if (new_mapped_indices[closest_vert_bindings[i]] == -1) { //Removed vertex
			is_loop = false;
			points_were_removed = true;
			continue;
		}
		new_bindings.push_back(new_mapped_indices[closest_vert_bindings[i]]);
		if (points_were_removed && no_tracked_point_yet) {
			first_included_after_remove = new_bindings.size() - 1;
			no_tracked_point_yet = false;
		}
	}

	if (new_bindings.size() == 0) {
		return false; //Stroke ceases to exist.
	}

	if (first_included_after_remove != -1) { //Check that we didn't remove ALL points in the stroke and that we actually need to reorder (because point 0 was removed)
		rotate(new_bindings.begin(), new_bindings.begin() + first_included_after_remove, new_bindings.end()); //Makes the stroke "continuous" again
	}

	new_bindings.push_back(new_bindings[0]); //Make the vertex bindings looped again
	set_closest_vert_bindings(new_bindings);

	return true;
}

void Stroke::undo_stroke_add(Eigen::VectorXi& edge_boundary_markers, Eigen::VectorXi& sharp_edge, Eigen::VectorXi& vertex_is_fixed) {
	Eigen::MatrixXi EV, FE, EF;
	igl::edge_topology(*V, *F, EV, FE, EF);

	int start, end, equal_pos;
	Eigen::VectorXi col1Equals, col2Equals;

	for (int i = 0; i < stroke_edges.rows() - !is_loop; i++) {
		start = stroke_edges(i, 0);
		end = stroke_edges(i, 1);
		col1Equals = EV.col(0).cwiseEqual(min(start, end)).cast<int>();
		col2Equals = EV.col(1).cwiseEqual(max(start, end)).cast<int>();
		int maxval = (col1Equals + col2Equals).maxCoeff(&equal_pos); //Find the row that contains both vertices of this edge
		if (maxval == 2) { //Make sure that the 2 sequential stroke3DPoints are actually connected (because we might have interrupted strokes)
			edge_boundary_markers[equal_pos] = 0;
			sharp_edge[equal_pos] = 0; //Unset sharp edges (e.g. edges that were previously sharp, now smooth, but they're not entirely removed)
		}
	}

	vector<vector<int>> neighbors;
	adjacency_list(*F, neighbors);
	for (int i = 0; i < closest_vert_bindings.size(); i++) {
		vertex_is_fixed[closest_vert_bindings[i]] = 0;
		for (int j = 0; j < neighbors[closest_vert_bindings[i]].size(); j++) {
			start = closest_vert_bindings[i];
			end = neighbors[closest_vert_bindings[i]][j];
			col1Equals = EV.col(0).cwiseEqual(min(start, end)).cast<int>();
			col2Equals = EV.col(1).cwiseEqual(max(start, end)).cast<int>();
			(col1Equals + col2Equals).maxCoeff(&equal_pos); //Find the row that contains both vertices of this edge
			if (edge_boundary_markers[equal_pos] > 0) { //If the vertex is part of another boundary edge, then it stays fixed
				vertex_is_fixed[closest_vert_bindings[i]] = 1;
				break;
			}
		}
	}
}

void Stroke::switch_stroke_edges_type(Eigen::VectorXi& sharp_edge) {
	Eigen::MatrixXi EV, FE, EF;
	igl::edge_topology(*V, *F, EV, FE, EF);

	int start, end, equal_pos;
	Eigen::VectorXi col1Equals, col2Equals;

	for (int i = 0; i < stroke_edges.rows() - !is_loop; i++) {
		start = stroke_edges(i, 0);
		end = stroke_edges(i, 1);
		col1Equals = EV.col(0).cwiseEqual(min(start, end)).cast<int>();
		col2Equals = EV.col(1).cwiseEqual(max(start, end)).cast<int>();
		int maxval = (col1Equals + col2Equals).maxCoeff(&equal_pos); //Find the row that contains both vertices of this edge
		if (maxval == 2) { //Make sure that the 2 sequential stroke3DPoints are actually connected (because we might have interrupted strokes)
			sharp_edge[equal_pos] = !sharp_edge[equal_pos];
		}
	}
}

//Only works for strokes that have valid stroke2DPoints (e.g. cut, extrusion base and extrusion silhouette strokes). Naive implementation in O(n^2)
bool Stroke::has_self_intersection(bool make_looped) {
	if (make_looped) {
		stroke2DPoints.conservativeResize(stroke2DPoints.rows() + 1, Eigen::NoChange);
		stroke2DPoints.bottomRows(1) = stroke2DPoints.row(0);
	}
	Eigen::RowVector2d p1, p2, p3, p4;
	for (int i = 0; i < stroke2DPoints.rows() - 1; i++) {
		p1 = stroke2DPoints.row(i);
		p2 = stroke2DPoints.row(i + 1);
		for (int j = i + 1; j < stroke2DPoints.rows() - 1; j++) {
			p3 = stroke2DPoints.row(j);
			p4 = stroke2DPoints.row(j + 1);
			if (line_segments_intersect(p1, p2, p3, p4)) {
				return true;
			}
		}
	}

	if (make_looped) {
		stroke2DPoints = stroke2DPoints.topRows(stroke2DPoints.rows() - 1);
	}
	return false;
}

bool Stroke::line_segments_intersect(Eigen::RowVector2d& p1, Eigen::RowVector2d& p2, Eigen::RowVector2d& p3, Eigen::RowVector2d& p4) {
	double a0, b0, c0, a1, b1, c1;
	a0 = p1[1] - p2[1];
	b0 = p2[0] - p1[0];
	c0 = p2[1] * p1[0] - p2[0] * p1[1];
	a1 = p3[1] - p4[1];
	b1 = p4[0] - p3[0];
	c1 = p4[1] * p3[0] - p4[0] * p3[1];

	if (((a0*p3[0] + b0 * p3[1] + c0)*(a0*p4[0] + b0 * p4[1] + c0) < 0) &&
		((a1*p1[0] + b1 * p1[1] + c1)*(a1*p2[0] + b1 * p2[1] + c1) < 0)) {
		return true;
	}
	else {
		return false;
	}
}

int Stroke::get_ID() {
	return stroke_ID;
}

Eigen::MatrixXd Stroke::get_V() const {
	return (*V);
}

Eigen::MatrixXi Stroke::get_F() const {
	return (*F);
}

Eigen::MatrixXd Stroke::get_stroke2DPoints() const {
	return stroke2DPoints;
}

Eigen::MatrixX3d Stroke::get3DPoints() {
	return stroke3DPoints;
}

Eigen::MatrixX3d Stroke::get3DPointsBack() {
	return stroke3DPointsBack;
}

Eigen::MatrixX3d Stroke::get_hand_pos() {
	return hand_pos_at_draw;
}

Eigen::MatrixXi Stroke::get_hit_faces() {
	return faces_hit;
}

Eigen::MatrixXi Stroke::get_stroke_edges() {
	return stroke_edges;
}

int Stroke::get_vertex_idx_for_point(int pt_idx) {
	return closest_vert_bindings[pt_idx];
}

vector<int> Stroke::get_closest_vert_bindings() {
	return closest_vert_bindings;
}

void Stroke::set3DPoints(Eigen::MatrixX3d new_3DPoints) {
	stroke3DPoints.resize(new_3DPoints.rows(), 3);
	stroke3DPoints = new_3DPoints;
}

void Stroke::set_closest_vert_bindings(vector<int> new_vert_bindings) {
	closest_vert_bindings.clear();
	closest_vert_bindings = new_vert_bindings;
}
