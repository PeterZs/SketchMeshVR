#ifdef _WIN32
#define NOMINMAX
#include <Windows.h>
#else
#include <unistd.h>
#include <curses.h>
#endif

#include <igl/readOFF.h>
#include <igl/edge_topology.h>
#include <igl/cat.h>
#include <igl/ray_mesh_intersect.h>
#include <iostream>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/boundary_loop.h>
#include <igl/harmonic.h>
#include <igl/map_vertices_to_circle.h>
#include <SketchMeshVR.h>
#include <Stroke.h>
#include "SurfaceSmoothing.h"
#include "MeshCut.h"
#include "CurveDeformation.h"
#include "MeshExtrusion.h"
#include <igl/opengl/OculusVR.h>


using namespace std;
using Viewer = igl::opengl::glfw::Viewer;
using OculusVR = igl::opengl::OculusVR;

Viewer viewer; 

// Vertex array, #V x3
Eigen::MatrixXd V(0,3);
// Face array, #F x3
Eigen::MatrixXi F(0,3);

//Per vertex indicator of whether vertex is on boundary (on boundary if == 1)
Eigen::VectorXi vertex_boundary_markers;
//Per vertex indicator of whether vertex is on original stroke (outline of shape) (on OG stroke if ==1)
Eigen::VectorXi part_of_original_stroke;
//Per edge indicator of whether the edge is sharp (if == 1 then sharp, otherwise smooth)
Eigen::VectorXi sharp_edge;
//Takes care of index mapping from before a cut/extrusion action to after (since some vertices are removed)
Eigen::VectorXi new_mapped_indices;

//General
enum ToolMode { DRAW, ADD, CUT, EXTRUDE, PULL, REMOVE, CHANGE, SMOOTH, NAVIGATE, TOGGLE, NONE, DEFAULT, FAIL }; //NONE is used to indicate that a button was released, whereas DEFAULT indicates that one of the toggle buttons was pressed within its cooldown period, FAIL is used to indicate that something went wrong (e.g. user clicked too far away for PULL)

ToolMode tool_mode = NAVIGATE;
ToolMode prev_tool_mode = NONE;
Stroke* initial_stroke;
Stroke* added_stroke;
Stroke* extrusion_base;
vector<Stroke> stroke_collection;

//For smoothing
int initial_smooth_iter = 30;

//For selecting vertices
int handleID = -1;


//Variables for pulling a curve (and removing added control curves)
int turnNr = 0;
bool dirty_boundary = false;
int closest_stroke_ID, prev_closest_stroke_ID;

//Keeps track of the stroke IDs
int next_added_stroke_ID = 2; //Start at 2 because marker 1 belongs to the original boundary


//Variables for removing a control curve
bool stroke_was_removed = false;
int remove_stroke_clicked = 0;

//Variables for cutting
bool cut_stroke_already_drawn = false;

//Variables for extrusion
bool extrusion_base_already_drawn = false;
SurfacePath base_surface_path;
Eigen::Matrix4f base_model, base_view, base_proj;
Eigen::Vector4f base_viewport;

bool button_A_is_set = false, button_B_is_set = false, button_thumb_is_set = false;

std::chrono::steady_clock::time_point _start_time, _end_time;
bool has_recentered = false;
bool draw_should_block = false;

void draw_all_strokes(){
	Eigen::MatrixXd added_points;
	viewer.data().set_points((Eigen::MatrixXd) initial_stroke->get3DPoints().block(0, 0, initial_stroke->get3DPoints().rows() - 1, 3), Eigen::RowVector3d(1, 0, 0)); //Display the original stroke points and clear all the rest. Don't take the last point
	viewer.data().set_edges(Eigen::MatrixXd(), Eigen::MatrixXi(), Eigen::RowVector3d(0, 0, 1)); //Clear the non-original stroke edges

	if (initial_stroke->is_loop) {
		viewer.data().set_stroke_points(igl::cat(1, (Eigen::MatrixXd) initial_stroke->get3DPoints().block(0, 0, initial_stroke->get3DPoints().rows() - 1, 3), (Eigen::MatrixXd) initial_stroke->get3DPoints().row(0))); //Create a loop and draw edges
	}
	else { //set_stroke_points always makes a loop, so don't use that when our stroke ain't a loop (anymore)
		added_points = initial_stroke->get3DPoints();
		viewer.data().add_edges(added_points.block(0, 0, added_points.rows() - 2, 3), added_points.block(1, 0, added_points.rows() - 2, 3), Eigen::RowVector3d(1, 0, 0));
	}

	int points_to_hold_back;
	for (int i = 0; i < stroke_collection.size(); i++) {
		added_points = stroke_collection[i].get3DPoints();
		points_to_hold_back = 1 + !stroke_collection[i].is_loop;
		viewer.data().add_points(added_points, stroke_collection[i].stroke_color);
		viewer.data().add_edges(added_points.block(0, 0, added_points.rows() - points_to_hold_back, 3), added_points.block(1, 0, added_points.rows() - points_to_hold_back, 3), stroke_collection[i].stroke_color);
	}
}

void draw_extrusion_base(){
	int points_to_hold_back;
	Eigen::MatrixXd added_points = extrusion_base->get3DPoints();
	points_to_hold_back = 1 + extrusion_base->is_loop;
	viewer.data().add_points(added_points, extrusion_base->stroke_color);
	viewer.data().add_edges(added_points.block(0, 0, added_points.rows() - points_to_hold_back, 3), added_points.block(1, 0, added_points.rows() - points_to_hold_back, 3), extrusion_base->stroke_color);
}

void select_dragging_handle(Eigen::Vector3f& pos) {
	double closest_dist = INFINITY;
	handleID = initial_stroke->selectClosestVertex(pos, closest_dist);
	double current_closest = closest_dist;
	closest_stroke_ID = -1;
	int tmp_handleID;
	for (int i = 0; i < stroke_collection.size(); i++) { //Additional strokes that cross the original stroke will never be selected as the pulled curve when the user clicks a vertex that also belongs to the original boundary, since their vertex positions are the same and we check for SMALLER distances
		tmp_handleID = stroke_collection[i].selectClosestVertex(pos, closest_dist); //Returns the index into stroke3DPoints of the closest point
		if ((closest_dist < current_closest) && (tmp_handleID != -1)) {
			current_closest = closest_dist;
			handleID = tmp_handleID;
			closest_stroke_ID = i;
		}
	}
}

ToolMode get_chosen_mode(OculusVR::ButtonCombo pressed) {
	if (pressed == OculusVR::ButtonCombo::GRIPTRIG) {
		if (button_B_is_set) {
			return PULL;
		} else {
			return DRAW;
		}
	}
	else if (pressed == OculusVR::ButtonCombo::GRIP) {
		if (button_A_is_set) {
			return EXTRUDE;
		} else {
			return CUT;
		}
	}
	else if (pressed == OculusVR::ButtonCombo::TRIG) {
		if (button_thumb_is_set) {
			return REMOVE;
		} else {
		//	return ADD;
			return SMOOTH;
		}
	}
	else if (pressed == OculusVR::ButtonCombo::A) {
		_start_time = _end_time; //Set timer 1 to be the previous time of when we came here
		_end_time = std::chrono::high_resolution_clock::now();
		auto timePast = std::chrono::duration_cast<std::chrono::nanoseconds>(_end_time - _start_time).count();
		if (timePast > 100000000) {
			button_A_is_set = !button_A_is_set;
			extrusion_base_already_drawn = false;
			return TOGGLE;
		}
		else {
			return DEFAULT;
		}
	}
	else if (pressed == OculusVR::ButtonCombo::B) {
		_start_time = _end_time; //Set timer 1 to be the previous time of when we came here
		_end_time = std::chrono::high_resolution_clock::now();
		auto timePast = std::chrono::duration_cast<std::chrono::nanoseconds>(_end_time - _start_time).count();
		if (timePast > 100000000) {
			button_B_is_set = !button_B_is_set;
			return TOGGLE;
		}
		else {
			return DEFAULT;
		}
	}
	else if (pressed == OculusVR::ButtonCombo::THUMB) {
		_start_time = _end_time; //Set timer 1 to be the previous time of when we came here
		_end_time = std::chrono::high_resolution_clock::now();
		auto timePast = std::chrono::duration_cast<std::chrono::nanoseconds>(_end_time - _start_time).count();
		if (timePast > 100000000) {
			button_thumb_is_set = !button_thumb_is_set;
			return TOGGLE;
		}
		else {
			return DEFAULT;
		}
	}
	else if (pressed == OculusVR::ButtonCombo::NONE) {
		return NONE;
	}
	else {
		return DEFAULT;
	}
}

void button_down(OculusVR::ButtonCombo pressed, Eigen::Vector3f& pos){
	ToolMode pressed_type = get_chosen_mode(pressed);

	if (pressed_type == DEFAULT) { //Default means that one of the toggle buttons was pressed again within the cooldown period. Do nothing
		return;
	}
	else if (pressed_type == TOGGLE) { //User was just switching between e.g. cut/extrude, don't do anything 
		prev_tool_mode = TOGGLE; //Needed for handling multithreading in VR_Viewer
		return;
	}
	else if (pressed_type == PULL || pressed_type == ADD || pressed_type == REMOVE || pressed_type == CUT || pressed_type == EXTRUDE) {
		if (initial_stroke->empty2D()) { //Don't go into these modes when there is no mesh yet
			prev_tool_mode = FAIL;
			return;
		}
	}
	else if (pressed_type == REMOVE) {
		if (stroke_collection.size() == 0) {
			return;
		}
	}

	tool_mode = pressed_type;

	Eigen::Vector3f pos_tmp = pos;
	Eigen::Vector3f last_eye_origin = viewer.oculusVR.get_last_eye_origin();
	pos_tmp[0] += last_eye_origin[0];
	pos_tmp[2] += last_eye_origin[2];
	if (tool_mode != DRAW && tool_mode != PULL) {
		Eigen::MatrixX3d LP(2, 3);
		LP.row(0) = pos_tmp.cast<double>();
		LP.row(1) = (pos_tmp + 1000 * viewer.oculusVR.get_right_touch_direction()).cast<double>();
		cout << "LP: " << LP.row(0) << endl;
		viewer.data().set_laser_points(LP);
		viewer.data().set_hand_point(pos_tmp.cast<double>().transpose(), Eigen::RowVector3d(0.5f, 0.5f, 0.5f));
	}
	else {
		viewer.data().set_laser_points(Eigen::MatrixXd(0, 0));
		viewer.data().set_hand_point(pos_tmp.cast<double>().transpose(), Eigen::RowVector3d(0.5f, 0.5f, 0.5f));
	}

	if (tool_mode == DRAW) { //Creating the first curve/mesh
		if (draw_should_block) { //User has been too close to first sample point (closing the stroke too much), so we're in blocked state till the buttons are released again
			return;
		}
		if (prev_tool_mode == NONE || prev_tool_mode == FAIL) {
			if (has_recentered) {
				viewer.oculusVR.set_start_action_view(viewer.core.get_view());
				stroke_collection.clear();
				next_added_stroke_ID = 2;
				initial_stroke->strokeReset();
				initial_stroke->addSegment(pos);
				prev_tool_mode = DRAW;
			}
			else {
				cut_stroke_already_drawn = false;
				extrusion_base_already_drawn = false;
				viewer.data().clear();
				viewer.oculusVR.request_recenter();
				has_recentered = true;
			}
		}
		else if (prev_tool_mode == DRAW) {
			//We had already started drawing, continue
			draw_should_block = initial_stroke->addSegment(pos);
			return;
		}
		else if (prev_tool_mode == CUT || prev_tool_mode == EXTRUDE || prev_tool_mode == ADD || prev_tool_mode == REMOVE) {
			// One of the modes was called while the user had only partially pressed the DRAW button combo. Stop drawing in Viewer and wait a round
			prev_tool_mode = NONE;
			viewer.update_screen_while_computing = false;
		}
	}
	else if (tool_mode == ADD) {
		if (prev_tool_mode == NONE || prev_tool_mode == FAIL) { //Adding a new control curve onto an existing mesh
			added_stroke = new Stroke(V, F, viewer, next_added_stroke_ID);
			next_added_stroke_ID++;
			added_stroke->addSegmentAdd(pos);
			prev_tool_mode = ADD;
		}
		else if (prev_tool_mode == ADD) {
			added_stroke->addSegmentAdd(pos);
		}
	}
	else if (tool_mode == REMOVE) {
		if (prev_tool_mode == NONE || prev_tool_mode == FAIL) {
			Eigen::Vector3f hit_pos;
			vector<igl::Hit> hits;
			Eigen::Vector3f pos_tmp = pos;
			Eigen::Vector3f last_eye_origin = viewer.oculusVR.get_last_eye_origin();
			pos_tmp[0] += last_eye_origin[0];
			pos_tmp[2] += last_eye_origin[2];

			if (igl::ray_mesh_intersect(pos_tmp, viewer.oculusVR.get_right_touch_direction(), V, F, hits)) { //Intersect the ray from the Touch controller with the mesh to get the 3D point
				hit_pos = (V.row(F(hits[0].id, 0))*(1.0 - hits[0].u - hits[0].v) + V.row(F(hits[0].id, 1))*hits[0].u + V.row(F(hits[0].id, 2))*hits[0].v).cast<float>();
			}
			else { //Hand ray did not intersect mesh
				prev_tool_mode = FAIL;
				return; 
			}

			double closest_dist = INFINITY;
			double current_closest = closest_dist;
			int tmp_handleID, closest_stroke_idx;
			handleID = -1;
			for (int i = 0; i < stroke_collection.size(); i++) {
				tmp_handleID = stroke_collection[i].selectClosestVertex(hit_pos, closest_dist);
				if ((closest_dist < current_closest) && (tmp_handleID != -1)) {
					current_closest = closest_dist;
					handleID = tmp_handleID;
					closest_stroke_ID = stroke_collection[i].get_ID();
					closest_stroke_idx = i;
				}
			}

			if (handleID == -1) {//User clicked too far from any of the stroke vertices
				prev_tool_mode = FAIL;
				return;
			}
			if (closest_stroke_ID == prev_closest_stroke_ID) {
				remove_stroke_clicked++;
			}
			else {
				remove_stroke_clicked = 1; //Start from 1
				prev_closest_stroke_ID = closest_stroke_ID;
			}

			//Redraw the original stroke and all added strokes, where the selected stroke is drawn in black.
			Eigen::MatrixXd init_points = initial_stroke->get3DPoints();
			viewer.data().set_points(init_points.topRows(init_points.rows() - 1), Eigen::RowVector3d(1, 0, 0));
			Eigen::MatrixXd added_points = stroke_collection[closest_stroke_idx].get3DPoints();
			viewer.data().add_points(added_points.topRows(added_points.rows() - 1), Eigen::RowVector3d(0, 0, 0));
			for (int i = 0; i < stroke_collection.size(); i++) {
				if (stroke_collection[i].get_ID() == closest_stroke_ID) {
					continue;
				}
				added_points = stroke_collection[i].get3DPoints();
				viewer.data().add_points(added_points.topRows(added_points.rows() - 1), stroke_collection[i].stroke_color);
			}

			if (remove_stroke_clicked == 2) { //Mechanism to force the user to click twice on the same stroke before removing it (safeguard)
				stroke_was_removed = true;
				stroke_collection[closest_stroke_idx].undo_stroke_add(vertex_boundary_markers); //Sets the vertex_boundary_markers for the vertices of this stroke to 0 again
				stroke_collection.erase(stroke_collection.begin() + closest_stroke_idx);
				remove_stroke_clicked = 0; //Reset
			}

			prev_tool_mode = REMOVE;
		}
		if (prev_tool_mode == REMOVE) {
			return; //For REMOVE we only take action upon button press and release 
		}
	}
	else if (tool_mode == PULL) { //Dragging an existing curve
		if (prev_tool_mode == NONE || prev_tool_mode == FAIL) { //Also allow to go to pull after ADD because sometimes the buttons are hard to differentiate
			Eigen::Vector3f last_eye_origin = viewer.oculusVR.get_last_eye_origin();
			pos[0] += last_eye_origin[0];
			pos[2] += last_eye_origin[2];
			select_dragging_handle(pos);

			if (handleID == -1) {//User clicked too far from any of the stroke vertices
				prev_tool_mode = FAIL;
				return;
			}
			if (closest_stroke_ID == -1) {
				CurveDeformation::startPullCurve(*initial_stroke, handleID);
			}
			else {
				CurveDeformation::startPullCurve(stroke_collection[closest_stroke_ID], handleID);
			}
			prev_tool_mode = PULL;
		}
		else if (prev_tool_mode == PULL) {
			Eigen::Vector3f last_eye_origin = viewer.oculusVR.get_last_eye_origin();
            pos[0] += last_eye_origin[0];
            pos[2] += last_eye_origin[2];
			if (turnNr == 0) { 
				CurveDeformation::pullCurve(pos.transpose().cast<double>(), V, part_of_original_stroke);
				SurfaceSmoothing::smooth(V, F, vertex_boundary_markers, part_of_original_stroke, new_mapped_indices, sharp_edge, dirty_boundary);

				turnNr++;

				initial_stroke->update_Positions(V);
				for (int i = 0; i < stroke_collection.size(); i++) {
					stroke_collection[i].update_Positions(V);
				}

				viewer.data().set_mesh(V, F);
				draw_all_strokes();
			}
			else {
				turnNr++;
				if (turnNr == 4) {//increase this number to smooth less often
					turnNr = 0;
				}
			}
		}
		else if (prev_tool_mode == CUT || prev_tool_mode == EXTRUDE || prev_tool_mode == ADD || prev_tool_mode == REMOVE) {
			// One of the modes was called while the user had only partially pressed the PULL button combo. Stop drawing in Viewer and wait a round
			prev_tool_mode = NONE;
			viewer.update_screen_while_computing = false;
		}
	}
	else if (tool_mode == CUT) {
		if (prev_tool_mode == NONE || prev_tool_mode == FAIL) {
			prev_tool_mode = CUT;
			if (cut_stroke_already_drawn) { //clicked while cut stroke already drawn
				return;
			}
			//clicked with no cut stroke drawn yet
			added_stroke = new Stroke(V, F, viewer, next_added_stroke_ID);
			viewer.oculusVR.set_start_action_view(viewer.core.get_view());
			next_added_stroke_ID++;
			added_stroke->addSegmentCut(pos);
		}
		else if (prev_tool_mode == CUT) {
			if (!cut_stroke_already_drawn) {
				added_stroke->addSegmentCut(pos);
			}
		}
	}
	else if (tool_mode == EXTRUDE) {
		if (prev_tool_mode == NONE || prev_tool_mode == FAIL) {
			prev_tool_mode = EXTRUDE;
			if (extrusion_base_already_drawn) { //clicked while the extrude base was already drawn
				added_stroke = new Stroke(V, F, viewer, next_added_stroke_ID);
				next_added_stroke_ID++;
				added_stroke->addSegmentExtrusionSilhouette(pos);
			}
			else { //clicked with no extrude base yet
				extrusion_base = new Stroke(V, F, viewer, next_added_stroke_ID);
				next_added_stroke_ID++;
				extrusion_base->addSegmentExtrusionBase(pos);
			}
		}
		else if (prev_tool_mode == EXTRUDE) {
			if (extrusion_base_already_drawn) {
				added_stroke->addSegmentExtrusionSilhouette(pos);
			}
			else {
				extrusion_base->addSegmentExtrusionBase(pos);
			}
		}
	}
	else if (tool_mode == SMOOTH) {
		prev_tool_mode = SMOOTH;
	}
	else if (tool_mode == NONE) {	//Have to finish up as if we're calling mouse_up()
		has_recentered = false;
		if (prev_tool_mode == NONE) {
			return;
		}
		else if (prev_tool_mode == DRAW) {
			prev_tool_mode = NONE;
			draw_should_block = false; 

			if (initial_stroke->toLoop()) {//Returns false if the stroke only consists of 1 point (user just clicked)

				initial_stroke->generate3DMeshFromStroke(vertex_boundary_markers, part_of_original_stroke, V, F);
				
				if (!igl::is_edge_manifold(F)) { //Check if the drawn stroke results in an edge-manifold mesh, otherwise sound a beep and revert
#ifdef _WIN32
					Beep(500, 200);
#else
					beep();
#endif				
					Eigen::MatrixXd drawn_points = initial_stroke->get3DPoints();
					initial_stroke->strokeReset();
					vertex_boundary_markers.resize(0);
					part_of_original_stroke.resize(0);
					dirty_boundary = true;

					viewer.data().add_edges(drawn_points.block(0, 0, drawn_points.rows() - 1, 3), drawn_points.block(1, 0, drawn_points.rows() - 1, 3), Eigen::RowVector3d(0, 0, 0)); //Display the stroke in black to show that it went wrong
					prev_tool_mode = NONE;
					viewer.update_screen_while_computing = false;
					return;
				}

				Eigen::MatrixXi EV, FE, EF;
				igl::edge_topology(V, F, EV, FE, EF);
				sharp_edge.resize(EV.rows());
				sharp_edge.setZero(); //Set all edges to smooth after initial draw

				dirty_boundary = true;
				for (int i = 0; i < initial_smooth_iter; i++) {
					SurfaceSmoothing::smooth(V, F, vertex_boundary_markers, part_of_original_stroke, new_mapped_indices, sharp_edge, dirty_boundary);
				}

				initial_stroke->update_Positions(V);

				viewer.data().clear();
				//viewer.data().set_face_based(true);
				viewer.data().set_mesh(V, F);

				//Overlay the drawn stroke
				int strokeSize = (vertex_boundary_markers.array() > 0).count();
				Eigen::MatrixXd strokePoints = V.block(0, 0, strokeSize, 3);
				viewer.data().set_points(strokePoints, Eigen::RowVector3d(1, 0, 0)); //Displays dots
				viewer.data().set_stroke_points(igl::cat(1, strokePoints, (Eigen::MatrixXd) V.row(0)));


			}
		}
		else if (prev_tool_mode == ADD) {
			dirty_boundary = true;
			if (!added_stroke->has_points_on_mesh) {
				viewer.update_screen_while_computing = false;
				return;
			}
			added_stroke->snap_to_vertices(vertex_boundary_markers);
			stroke_collection.push_back(*added_stroke);
			draw_all_strokes();
		}
		else if (prev_tool_mode == REMOVE && stroke_was_removed) { //Only redraw if we actually removed a stroke (otherwise we draw unnecessary)
			stroke_was_removed = false; //Reset
			dirty_boundary = true;

			draw_all_strokes();
		}
		else if (prev_tool_mode == PULL && handleID != -1) {
			for (int i = 0; i < 6; i++) {
				SurfaceSmoothing::smooth(V, F, vertex_boundary_markers, part_of_original_stroke, new_mapped_indices, sharp_edge, dirty_boundary);
			}

			for (int i = 0; i < stroke_collection.size(); i++) {
				stroke_collection[i].update_Positions(V);
			}
			viewer.data().clear();
			viewer.data().set_mesh(V, F);
			draw_all_strokes();
		}
		else if (prev_tool_mode == CUT) {
			if (!added_stroke->has_points_on_mesh) {
				viewer.update_screen_while_computing = false;
				return;
			}
			if (cut_stroke_already_drawn) { //User had already drawn the cut stroke and has now clicked/drawn the final stroke for removing the part
				dirty_boundary = true;

				added_stroke->prepend_first_point();
				added_stroke->append_final_point();
				added_stroke->toLoop();
				
				bool cut_success = MeshCut::cut(V, F, vertex_boundary_markers, part_of_original_stroke, new_mapped_indices, sharp_edge, *added_stroke);
				if (!cut_success) { //Catches the cases that the cut removes all mesh vertices/faces and when the first/last cut point aren't correct (face -1 in SurfacePath)
					cut_stroke_already_drawn = false;
					next_added_stroke_ID--; //Undo ID increment since stroke didn't actually get pushed back
					prev_tool_mode = NONE;
#ifdef _WIN32
					Beep(500, 200);
#else
					beep();
#endif	

					draw_all_strokes(); //Will remove the pink cut stroke
					Eigen::MatrixXd drawn_points = added_stroke->get3DPoints();
					viewer.data().add_edges(drawn_points.block(0, 0, drawn_points.rows() - 1, 3), drawn_points.block(1, 0, drawn_points.rows() - 1, 3), Eigen::RowVector3d(0, 0, 0)); //Display the stroke in black to show that it went wrong
					//draw_all_strokes();
					viewer.update_screen_while_computing = false;
					return;
				}
				stroke_collection.push_back(*added_stroke);

				initial_stroke->update_vert_bindings(new_mapped_indices, vertex_boundary_markers); //Don't test if the initial one dies

				int nr_removed = 0, original_collection_size = stroke_collection.size();
				for (int i = 0; i < original_collection_size; i++) {
					if (!stroke_collection[i - nr_removed].update_vert_bindings(new_mapped_indices, vertex_boundary_markers)) {
						//Stroke dies, don't need to do stroke.undo_stroke_add, cause all its vertices also cease to exist
						stroke_collection.erase(stroke_collection.begin() + i - nr_removed);
						nr_removed++;
						dirty_boundary = true;
						continue; //Go to the next stroke, don't update this ones' positions
					}
				}

				for (int i = 0; i < 10; i++) {
					SurfaceSmoothing::smooth(V, F, vertex_boundary_markers, part_of_original_stroke, new_mapped_indices, sharp_edge, dirty_boundary);
				}

				//Update the stroke positions after smoothing, in case their positions have changed (although they really shouldn't)
				initial_stroke->update_Positions(V);
				for (int i = 0; i < stroke_collection.size(); i++) {
					stroke_collection[i].update_Positions(V);
				}

				viewer.data().clear();
				viewer.data().set_mesh(V, F);

				cut_stroke_already_drawn = false; //Reset
				draw_all_strokes();
			}
			else { //We're finished drawing the cut stroke, prepare for when user draws the final stroke to remove the part
				cut_stroke_already_drawn = true;
			}
		}
		else if (prev_tool_mode == EXTRUDE) {
			if (extrusion_base_already_drawn) { //User has drawn the silhouette stroke for extrusion
				dirty_boundary = true;
				added_stroke->toLoop();
				added_stroke->is_loop = false; //Set to false manually, because we don't want curveDeformation to consider it as a loop (but we do need looped 3DPoints)

				bool success_extrude = MeshExtrusion::extrude_main(V, F, vertex_boundary_markers, part_of_original_stroke, new_mapped_indices, sharp_edge, base_surface_path, *added_stroke, *extrusion_base, base_model, base_view, base_proj, base_viewport);
				if (!success_extrude) { //Catches the case that the extrusion base removes all faces/vertices
					prev_tool_mode = NONE;
					next_added_stroke_ID -= 2;
					extrusion_base_already_drawn = false;
#ifdef _WIN32
					Beep(700, 200);
#else
					beep();
#endif	
					draw_all_strokes(); //Will remove the drawn base & silhouette strokes
					Eigen::MatrixXd drawn_points = extrusion_base->get3DPoints();
					viewer.data().add_edges(drawn_points.block(0, 0, drawn_points.rows() - 1, 3), drawn_points.block(1, 0, drawn_points.rows() - 1, 3), Eigen::RowVector3d(0, 0, 0)); //Display the stroke in black to show that it went wrong
					drawn_points = added_stroke->get3DPoints();
					viewer.data().add_edges(drawn_points.block(0, 0, drawn_points.rows() - 2, 3), drawn_points.block(1, 0, drawn_points.rows() - 2, 3), Eigen::RowVector3d(0, 0, 0)); //Display the stroke in black to show that it went wrong
					viewer.update_screen_while_computing = false;
					return;
				}


				stroke_collection.push_back(*extrusion_base);
				stroke_collection.push_back(*added_stroke);

				initial_stroke->update_vert_bindings(new_mapped_indices, vertex_boundary_markers);

				int nr_removed = 0, original_collection_size = stroke_collection.size();
				for (int i = 0; i < original_collection_size - 1; i++) { //Skip the newly added stroke since it is already updated inside extrude_main()
					if (!stroke_collection[i - nr_removed].update_vert_bindings(new_mapped_indices, vertex_boundary_markers)) {
						//Stroke dies, don't need to do stroke.undo_stroke_add, cause all its vertices also cease to exist (or in case of non-loop strokes that have a middle portion removed, the undo_stroke_add is done inside of the update_vert_bindings
						stroke_collection.erase(stroke_collection.begin() + i - nr_removed);
						nr_removed++;
						dirty_boundary = true;
						continue; //Go to the next stroke, don't update this ones' positions
					}
				}


				for (int i = 0; i < 10; i++) {
					SurfaceSmoothing::smooth(V, F, vertex_boundary_markers, part_of_original_stroke, new_mapped_indices, sharp_edge, dirty_boundary);
				}

				//Update the stroke positions after smoothing, in case their positions have changed (although they really shouldn't)
				initial_stroke->update_Positions(V);
				for (int i = 0; i < stroke_collection.size() - 1; i++) {
					stroke_collection[i].update_Positions(V);
				}


				viewer.data().clear();
				viewer.data().set_mesh(V, F);

				draw_all_strokes();

				extrusion_base_already_drawn = false;
			}
			else { //mouse released after extrusion base drawn
				if (!extrusion_base->has_points_on_mesh) {
					viewer.update_screen_while_computing = false;
					prev_tool_mode = NONE;
					return;
				}

				dirty_boundary = true;
				extrusion_base->toLoop();

				bool succes_extrude_prepare = MeshExtrusion::extrude_prepare(*extrusion_base, base_surface_path); //Don't need to update all strokes here, since it didn't remove any vertices
				if (!succes_extrude_prepare) { //Catches the case that face == -1 in SurfacePath
#ifdef _WIN32
					Beep(900, 200);
#else
					beep();
#endif	
					next_added_stroke_ID--;
					prev_tool_mode = NONE;
					draw_all_strokes(); //Removes the drawn base stroke
					Eigen::MatrixXd drawn_points = extrusion_base->get3DPoints();
					viewer.data().add_edges(drawn_points.block(0, 0, drawn_points.rows() - 1, 3), drawn_points.block(1, 0, drawn_points.rows() - 1, 3), Eigen::RowVector3d(0, 0, 0)); //Display the stroke in black to show that it went wrong
					viewer.update_screen_while_computing = false;
					extrusion_base_already_drawn = false;
					return;
				}

				extrusion_base_already_drawn = true;

				base_model = extrusion_base->viewer.core.get_model();
				base_view = extrusion_base->viewer.core.get_view();
				base_proj = extrusion_base->viewer.core.get_proj();
				base_viewport = extrusion_base->viewer.core.viewport;

				draw_all_strokes();
				draw_extrusion_base(); //Need to draw the extrusion base separately, since it isn't added to the stroke_collection yet.

			}
		}
		else if (prev_tool_mode == TOGGLE || prev_tool_mode == FAIL) {
			viewer.update_screen_while_computing = false;
		}
		else if (prev_tool_mode == SMOOTH) {
			SurfaceSmoothing::smooth(V, F, vertex_boundary_markers, part_of_original_stroke, new_mapped_indices, sharp_edge, dirty_boundary);
			for (int i = 0; i < stroke_collection.size(); i++) {
				stroke_collection[i].update_Positions(V);
			}
			viewer.data().clear();
			viewer.data().set_mesh(V, F);
			draw_all_strokes();
		}


		prev_tool_mode = NONE;
		viewer.update_screen_while_computing = false;
		return;
	}

	return;
}


bool callback_load_mesh(Viewer& viewer, string filename, Eigen::MatrixXd& V_floor, Eigen::MatrixXi& F_floor) {
		igl::readOFF(filename, V, F);
		viewer.data().clear();
		viewer.data().set_mesh(V_floor, F_floor);
		viewer.append_mesh();
		viewer.data().set_mesh(V, F);

		std::cout << filename.substr(filename.find_last_of("/") + 1) << endl;
	return true;
}

int main(int argc, char *argv[]) {
	//Init stroke selector
	initial_stroke = new Stroke(V, F, viewer, 0);

	Eigen::MatrixXd V_floor(4, 3);
	V_floor.row(0) << -100, 0, -100;
	V_floor.row(1) << 100, 0, -100;
	V_floor.row(2) << 100, 0, 100;
	V_floor.row(3) << -100, 0, 100;

	Eigen::MatrixXi F_floor(2, 3);
	F_floor.row(0) << 0, 3, 1;
	F_floor.row(1) << 3, 2, 1;

	if (argc == 2) {
		// Read mesh
		igl::readOFF(argv[1], V, F);
		callback_load_mesh(viewer, argv[1], V_floor, F_floor);
	}
	else {
		
		viewer.data().set_mesh(V_floor, F_floor);
		Eigen::MatrixXd V_uv(V_floor.rows(), 2);
		V_uv.col(0) = V_floor.col(0);
		V_uv.col(1) = V_floor.col(2);
		V_uv.col(0) = V_uv.col(0).array() - V_uv.col(0).minCoeff();
		V_uv.col(0) = V_uv.col(0).array() / V_uv.col(0).maxCoeff();
		V_uv.col(1) = V_uv.col(1).array() - V_uv.col(1).minCoeff();
		V_uv.col(1) = V_uv.col(1).array() / V_uv.col(1).maxCoeff();
		V_uv = V_uv.array() * 40;

		viewer.data().set_uv(V_uv); 
		viewer.data().show_texture = false; //TODO: texture turned off for now. Due to problems with anti-aliasing
		viewer.append_mesh();
		viewer.data().set_mesh(V, F);
		//viewer.load_mesh_from_file("../data/cube.off");
	}

	CurveDeformation::smooth_deform_mode = true;
	viewer.init_oculus();
	viewer.oculusVR.callback_button_down = button_down;
	viewer.data().point_size = 15;
	viewer.launch_oculus();
}

