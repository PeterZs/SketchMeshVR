#ifdef _WIN32
#define NOMINMAX
#include <Windows.h>
#else
#include <unistd.h>
#endif

#include <igl/readOFF.h>
#include <igl/read_triangle_mesh.h>
#include <igl/edge_topology.h>
#include <igl/cat.h>
#include <igl/ray_mesh_intersect.h>
#include <iostream>
#include <unordered_map>
#include <igl/viewer/VR_Viewer.h>
#include <SketchMeshVR.h>
#include <Stroke.h>
#include "SurfaceSmoothing.h"
#include "MeshCut.h"
#include "CurveDeformation.h"
#include "MeshExtrusion.h"


using namespace std;
using ViewerVR = igl::viewer::VR_Viewer;

ViewerVR viewervr; //TODO: I made this a global variable so we don't need to send a viewerVR into the callback for button down, which was giving problems due to copy constructors of classes involving atomic variables. Document any possible side effects of this (not sure if there was any benefit of sending in the viewerVR anyway)


// Vertex array, #V x3
Eigen::MatrixXd V(0,3);
// Face array, #F x3
Eigen::MatrixXi F(0,3);
// Per face normals, #F x3
Eigen::MatrixXd N_Faces;

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

//Mouse interaction
bool skip_standardcallback = false;
int down_mouse_x = -1, down_mouse_y = -1;
bool mouse_is_down = false; //We need this due to mouse_down not working in the nanogui menu, whilst mouse_up does work there
bool hand_has_moved = false;
Eigen::Vector3f prev_pos = Eigen::Vector3f::Zero();

//For smoothing
int initial_smooth_iter = 8;

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

bool button_A_is_set = false;
bool button_B_is_set = false;
bool button_thumb_is_set = false;

std::chrono::steady_clock::time_point _start_time, _end_time;
bool has_recentered = false;

void draw_all_strokes(){
	Eigen::MatrixXd added_points;
	viewervr.data.set_points((Eigen::MatrixXd) initial_stroke->get3DPoints().block(0, 0, initial_stroke->get3DPoints().rows() - 1, 3), Eigen::RowVector3d(1, 0, 0)); //Display the original stroke points and clear all the rest. Don't take the last point
	viewervr.data.set_edges(Eigen::MatrixXd(), Eigen::MatrixXi(), Eigen::RowVector3d(0, 0, 1)); //Clear the non-original stroke edges

	if (initial_stroke->is_loop) {
		viewervr.data.set_stroke_points(igl::cat(1, (Eigen::MatrixXd) initial_stroke->get3DPoints().block(0, 0, initial_stroke->get3DPoints().rows() - 1, 3), (Eigen::MatrixXd) initial_stroke->get3DPoints().row(0))); //Create a loop and draw edges
	}
	else { //set_stroke_points always makes a loop, so don't use that when our stroke ain't a loop (anymore)
		added_points = initial_stroke->get3DPoints();
		viewervr.data.add_edges(added_points.block(0, 0, added_points.rows() - 2, 3), added_points.block(1, 0, added_points.rows() - 2, 3), Eigen::RowVector3d(1, 0, 0));
	}

	int points_to_hold_back;
	for (int i = 0; i < stroke_collection.size(); i++) {
		added_points = stroke_collection[i].get3DPoints();
		points_to_hold_back = 1 + !stroke_collection[i].is_loop;
		viewervr.data.add_points(added_points, stroke_collection[i].stroke_color);
		viewervr.data.add_edges(added_points.block(0, 0, added_points.rows() - points_to_hold_back, 3), added_points.block(1, 0, added_points.rows() - points_to_hold_back, 3), stroke_collection[i].stroke_color);
	}
}

void draw_extrusion_base(){
	int points_to_hold_back;
	Eigen::MatrixXd added_points = extrusion_base->get3DPoints();
	points_to_hold_back = 1 + extrusion_base->is_loop;
	viewervr.data.add_points(added_points, extrusion_base->stroke_color);
	viewervr.data.add_edges(added_points.block(0, 0, added_points.rows() - points_to_hold_back, 3), added_points.block(1, 0, added_points.rows() - points_to_hold_back, 3), extrusion_base->stroke_color);
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

ToolMode get_chosen_mode(ViewerVR::ButtonCombo pressed) {
	if (pressed == ViewerVR::ButtonCombo::GRIPTRIG) {
		if (button_B_is_set) {
			cout << "pull mode" << endl;
			return PULL;
		} else {
			cout << "draw" << endl;
			return DRAW;
		}
	}
	else if (pressed == ViewerVR::ButtonCombo::GRIP) {
		if (button_A_is_set) {
			cout << "extrud" << endl;
			return EXTRUDE;
		} else {
			cout << "cut mode" << endl;

			return CUT;
		}
	}
	else if (pressed == ViewerVR::ButtonCombo::TRIG) {
		if (button_thumb_is_set) {
			cout << "remove mode" << endl;
			return REMOVE;
		} else {
			cout << "Add mode" << endl;
			return ADD;
		}
	}
	else if (pressed == ViewerVR::ButtonCombo::A) {
		_start_time = _end_time; //Set timer 1 to be the previous time of when we came here
		_end_time = std::chrono::high_resolution_clock::now();
		auto timePast = std::chrono::duration_cast<std::chrono::nanoseconds>(_end_time - _start_time).count();
		if (timePast > 100000000) {
			button_A_is_set = !button_A_is_set;
			cout << "Switching A " << endl;
			return TOGGLE;
		}
	}
	else if (pressed == ViewerVR::ButtonCombo::B) {
		_start_time = _end_time; //Set timer 1 to be the previous time of when we came here
		_end_time = std::chrono::high_resolution_clock::now();
		auto timePast = std::chrono::duration_cast<std::chrono::nanoseconds>(_end_time - _start_time).count();
		if (timePast > 100000000) {
			button_B_is_set = !button_B_is_set;
			cout << "switching buttonB" << endl;
			return TOGGLE;
		}
	}
	else if (pressed == ViewerVR::ButtonCombo::THUMB) {
		_start_time = _end_time; //Set timer 1 to be the previous time of when we came here
		_end_time = std::chrono::high_resolution_clock::now();
		auto timePast = std::chrono::duration_cast<std::chrono::nanoseconds>(_end_time - _start_time).count();
		if (timePast > 100000000) {
			button_thumb_is_set = !button_thumb_is_set;
			cout << "switching thumb button" << endl;
			return TOGGLE;
		}
	}
	else if (pressed == ViewerVR::ButtonCombo::NONE) {
		cout << "none" << endl;
		return NONE;
	}
	else {
		cout << "default" << endl;
		return DEFAULT;
	}
}

bool button_down(ViewerVR::ButtonCombo pressed, Eigen::Vector3f& pos){
	ToolMode pressed_type = get_chosen_mode(pressed);
	cout << "pressed type " << pressed_type << endl;
	if (pressed_type == DEFAULT) { //Default means that one of the toggle buttons was pressed again within the cooldown period. Do nothing
		return true;
	}
	if (pressed_type == TOGGLE) { //User was just switching between e.g. cut/extrude, don't do anything 
		prev_tool_mode = TOGGLE; //Needed for handling multithreading in VR_Viewer
		return true;
	}

/*	if (!((pos - prev_pos).isZero())) {
		hand_has_moved = true;
	}
	prev_pos = pos;*/

	if (pressed_type == PULL || pressed_type == ADD || pressed_type == CUT || pressed_type == EXTRUDE) {
		if (initial_stroke->empty2D()) { //Don't go into these modes when there is no mesh yet
			cout << "exiting due to initial stroke's empty2d" << endl;
			return true;
		}
	}
	if (pressed_type == REMOVE) {
		if (stroke_collection.size() == 0) {
			return true;
		}
		
		//	remove_stroke_clicked = 0; //Reset because we might be left with a single click from the last round. 
		//For now commented out because this gets called every iteration, whereas for a mouse it only got called when the button was pressed/released. Would need to add a var that tracks the last not-NONE mode
		
	}

	tool_mode = pressed_type;

	if (tool_mode != DRAW && tool_mode != PULL) {
		Eigen::MatrixX3d LP(2, 3);
		Eigen::Vector3f pos_tmp = pos;
		cout << "want to get eye pos" << endl;
		pos_tmp[0] += viewervr.get_current_eye_pos()[0];
		pos_tmp[2] += viewervr.get_current_eye_pos()[2];
		cout << "got eye pos" << endl;
		LP.row(0) = pos_tmp.cast<double>();
		LP.row(1) = (pos + 1000 * viewervr.get_right_touch_direction()).cast<double>();
		cout << "got touch dir" << endl;
		viewervr.data.set_laser_points(LP);
		viewervr.data.set_hand_point(pos_tmp.cast<double>().transpose(), Eigen::RowVector3d(0.5f, 0.5f, 0.5f));
		cout << "set overlay stuff" << endl;
	}
	else {
		Eigen::RowVector3f pos_tmp = pos;
		cout << "want pos for draw or pull" << endl;

		pos_tmp[0] += viewervr.get_current_eye_pos()[0];
		pos_tmp[2] += viewervr.get_current_eye_pos()[2];
		cout << "got pos for draw or pull" << endl;
		viewervr.data.set_hand_point(pos_tmp.cast<double>(), Eigen::RowVector3d(0.5f, 0.5f, 0.5f));
		viewervr.data.set_laser_points(Eigen::MatrixXd(0, 0));
		cout << "set overlay for draw or pull" << endl;
	}
	cout << "prev toolmode " << prev_tool_mode << endl;
	if (tool_mode == DRAW) { //Creating the first curve/mesh
		if (prev_tool_mode == NONE || prev_tool_mode == FAIL) {
			if (has_recentered) {
				viewervr.set_start_action_view(viewervr.corevr.get_view());
				stroke_collection.clear();
				next_added_stroke_ID = 2;
				initial_stroke->strokeReset();
				initial_stroke->strokeAddSegment(pos);
				prev_tool_mode = DRAW;
				skip_standardcallback = true;
			}
			else {
				viewervr.data.clear_without_floor();
				viewervr.request_recenter();
				has_recentered = true;
			}
		}
		else if (prev_tool_mode == DRAW) {
			//We had already started drawing, continue
			initial_stroke->strokeAddSegment(pos);
			return true;
		}
	}
	else if (tool_mode == ADD) {
		if (prev_tool_mode == NONE || prev_tool_mode == FAIL) { //Adding a new control curve onto an existing mesh
			added_stroke = new Stroke(V, F, viewervr, next_added_stroke_ID);
			next_added_stroke_ID++;
			added_stroke->strokeAddSegmentAdd(pos); //If the user starts outside of the mesh, consider the movement as navigation
			prev_tool_mode = ADD;
			skip_standardcallback = true;
		}
		else if (prev_tool_mode == ADD) {
			added_stroke->strokeAddSegmentAdd(pos);
			return true;
		}
	}
	else if (tool_mode == REMOVE) {
		if (prev_tool_mode == NONE || prev_tool_mode == FAIL) {
			Eigen::Vector3f hit_pos;
			vector<igl::Hit> hits;
			Eigen::Vector3f pos_tmp = pos;
			pos_tmp[0] += viewervr.get_current_eye_pos()[0];
			pos_tmp[2] += viewervr.get_current_eye_pos()[2];

			igl::ray_mesh_intersect(pos_tmp, viewervr.get_right_touch_direction(), V, F, hits); //Intersect the ray from the Touch controller with the mesh to get the 3D point
			hit_pos = (V.row(F(hits[0].id, 0))*(1.0 - hits[0].u - hits[0].v) + V.row(F(hits[0].id, 1))*hits[0].u + V.row(F(hits[0].id, 2))*hits[0].v).cast<float>();

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
				return false;
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
			viewervr.data.set_points(init_points.topRows(init_points.rows() - 1), Eigen::RowVector3d(1, 0, 0));
			Eigen::MatrixXd added_points = stroke_collection[closest_stroke_idx].get3DPoints();
			viewervr.data.add_points(added_points.topRows(added_points.rows() - 1), Eigen::RowVector3d(0, 0, 0));
			for (int i = 0; i < stroke_collection.size(); i++) {
				if (stroke_collection[i].get_ID() == closest_stroke_ID) {
					continue;
				}
				added_points = stroke_collection[i].get3DPoints();
				viewervr.data.add_points(added_points.topRows(added_points.rows() - 1), stroke_collection[i].stroke_color);
			}

			if (remove_stroke_clicked == 2) { //Mechanism to force the user to click twice on the same stroke before removing it (safeguard)
				stroke_was_removed = true;
				stroke_collection[closest_stroke_idx].undo_stroke_add(vertex_boundary_markers); //Sets the vertex_boundary_markers for the vertices of this stroke to 0 again
				stroke_collection.erase(stroke_collection.begin() + closest_stroke_idx);
				remove_stroke_clicked = 0; //Reset
			}

			prev_tool_mode = REMOVE;
			skip_standardcallback = true;
		}
		if (prev_tool_mode == REMOVE) {
			return true; //For REMOVE we only take action upon button press and release 
		}
	}
	else if (tool_mode == PULL) { //Dragging an existing curve
		if (prev_tool_mode == NONE || prev_tool_mode == ADD || prev_tool_mode == FAIL) { //Also allow to go to pull after ADD because sometimes the buttons are hard to differentiate
			pos[0] += viewervr.get_current_eye_pos()[0];
			pos[2] += viewervr.get_current_eye_pos()[2];
			select_dragging_handle(pos);

			if (handleID == -1) {//User clicked too far from any of the stroke vertices
				prev_tool_mode = FAIL;
				return false;
			}
			if (closest_stroke_ID == -1) {
				CurveDeformation::startPullCurve(*initial_stroke, handleID);
			}
			else {
				CurveDeformation::startPullCurve(stroke_collection[closest_stroke_ID], handleID);
			}
			prev_tool_mode = PULL;
			skip_standardcallback = true;
		}
		else if (prev_tool_mode == PULL) {
            pos[0] += viewervr.get_current_eye_pos()[0];
            pos[2] += viewervr.get_current_eye_pos()[2];
			if (turnNr == 0) { //increase the number to smooth less often
				CurveDeformation::pullCurve(pos.transpose().cast<double>(), V, part_of_original_stroke);

				SurfaceSmoothing::smooth(V, F, vertex_boundary_markers, part_of_original_stroke, new_mapped_indices, sharp_edge, dirty_boundary);

				turnNr++;

				initial_stroke->update_Positions(V);
				for (int i = 0; i < stroke_collection.size(); i++) {
					stroke_collection[i].update_Positions(V);
				}

				viewervr.data.set_mesh_with_floor(V, F);
				draw_all_strokes();

			}
			else {
				turnNr++;
				if (turnNr == 4) {
					turnNr = 0;
				}
			}

			
			return true;
		}

	}
	else if (tool_mode == CUT) {
		if (prev_tool_mode == NONE || prev_tool_mode == FAIL) {
			prev_tool_mode = CUT;
			if (cut_stroke_already_drawn) { //clicked while cut stroke already drawn
				cout << "return after second click" << endl;
				return true;
			}
			//clicked with no cut stroke drawn yet
			added_stroke = new Stroke(V, F, viewervr, next_added_stroke_ID);
			viewervr.set_start_action_view(viewervr.corevr.get_view());
			next_added_stroke_ID++;
			added_stroke->strokeAddSegmentCut(pos);
			skip_standardcallback = true;
		}
		else if (prev_tool_mode == CUT) {
			if (!cut_stroke_already_drawn) {
				added_stroke->strokeAddSegmentCut(pos);
			}
			return true;
		}
	}
	else if (tool_mode == EXTRUDE) {
		if (prev_tool_mode == NONE || prev_tool_mode == FAIL) {
			prev_tool_mode = EXTRUDE;
			if (extrusion_base_already_drawn) { //clicked while the extrude base was already drawn
				added_stroke = new Stroke(V, F, viewervr, next_added_stroke_ID);
				next_added_stroke_ID++;
				added_stroke->strokeAddSegmentExtrusionSilhouette(pos);
			}
			else { //clicked with no extrude base yet
				extrusion_base = new Stroke(V, F, viewervr, next_added_stroke_ID);
				next_added_stroke_ID++;
				extrusion_base->strokeAddSegmentExtrusionBase(pos);
			}
			skip_standardcallback = true;
		}
		else if (prev_tool_mode == EXTRUDE) {
			if (extrusion_base_already_drawn) {
				added_stroke->strokeAddSegmentExtrusionSilhouette(pos);
			}
			else {
				extrusion_base->strokeAddSegmentExtrusionBase(pos);
			}
			return true;
		}
	}
	else if (tool_mode == NONE) {	//Have to finish up as if we're calling mouse_up()
		has_recentered = false;
		if (prev_tool_mode == NONE) {
			return true;
		}

		else if (prev_tool_mode == DRAW) {
			prev_tool_mode = NONE;
			initial_stroke->strokeAddSegment(pos);

			if (initial_stroke->toLoop()) {//Returns false if the stroke only consists of 1 point (user just clicked)

				initial_stroke->generate3DMeshFromStroke(vertex_boundary_markers, part_of_original_stroke, V, F);
				
				Eigen::MatrixXi EV, FE, EF;
				igl::edge_topology(V, F, EV, FE, EF);
				sharp_edge.resize(EV.rows());
				sharp_edge.setZero(); //Set all edges to smooth after initial draw

				dirty_boundary = true;

				for (int i = 0; i < initial_smooth_iter; i++) {
					SurfaceSmoothing::smooth(V, F, vertex_boundary_markers, part_of_original_stroke, new_mapped_indices, sharp_edge, dirty_boundary);
				}

				initial_stroke->update_Positions(V);

				viewervr.data.clear_all();
				viewervr.data.set_face_based(true);
				viewervr.data.set_mesh_with_floor(V, F);

				//Overlay the drawn stroke
				int strokeSize = (vertex_boundary_markers.array() > 0).count();
				Eigen::MatrixXd strokePoints = V.block(0, 0, strokeSize, 3);
				viewervr.data.set_points(strokePoints, Eigen::RowVector3d(1, 0, 0)); //Displays dots
				viewervr.data.set_stroke_points(igl::cat(1, strokePoints, (Eigen::MatrixXd) V.row(0)));


			}
			hand_has_moved = false;
			skip_standardcallback = false;
		}
		else if (prev_tool_mode == ADD) {
			dirty_boundary = true;
			if (!added_stroke->has_points_on_mesh) {
				hand_has_moved = false;
				viewervr.draw_while_computing = false;
				return true;
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
		else if (prev_tool_mode == PULL && handleID != -1) { //TODO: do we need hand_has_moved logic?
			cout << "start finishing pull" << endl;
			for (int i = 0; i < 2; i++) {
				SurfaceSmoothing::smooth(V, F, vertex_boundary_markers, part_of_original_stroke, new_mapped_indices, sharp_edge, dirty_boundary);
			}

			for (int i = 0; i < stroke_collection.size(); i++) {
				stroke_collection[i].update_Positions(V);
			}
			cout << "finish pull before clear and set" << endl;
			viewervr.data.clear_all();
			viewervr.data.set_mesh_with_floor(V, F);
			cout << "finish pull before draw_all" << endl;
			draw_all_strokes();
			cout << "Finish pull after draw_all" << endl;
		}
		else if (prev_tool_mode == CUT) {
			if (!added_stroke->has_points_on_mesh) {
				hand_has_moved = false;
				viewervr.draw_while_computing = false;
				return true;
			}
			if (cut_stroke_already_drawn) { //User had already drawn the cut stroke and has now drawn the final stroke for removing the part
				dirty_boundary = true;

				added_stroke->prepend_first_point();
				added_stroke->append_final_point();
				added_stroke->toLoop();

				cout << added_stroke->get3DPoints() << endl << endl << added_stroke->get_stroke2DPoints() << endl << endl;
				MeshCut::cut(V, F, vertex_boundary_markers, part_of_original_stroke, new_mapped_indices, sharp_edge, *added_stroke);
				stroke_collection.push_back(*added_stroke);

				initial_stroke->update_vert_bindings(new_mapped_indices, vertex_boundary_markers);//Don't test if the initial one dies, cause then we have mayhem anyway? TODO

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

				for (int i = 0; i < 2; i++) {
					SurfaceSmoothing::smooth(V, F, vertex_boundary_markers, part_of_original_stroke, new_mapped_indices, sharp_edge, dirty_boundary);
				}

				//Update the stroke positions after smoothing, in case their positions have changed (although they really shouldn't)
				initial_stroke->update_Positions(V);
				for (int i = 0; i < stroke_collection.size(); i++) {
					stroke_collection[i].update_Positions(V);
				}

				viewervr.data.clear_all();
				viewervr.data.set_mesh_with_floor(V, F);

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
				MeshExtrusion::extrude_main(V, F, vertex_boundary_markers, part_of_original_stroke, new_mapped_indices, sharp_edge, base_surface_path, *added_stroke, *extrusion_base, base_model, base_view, base_proj, base_viewport);
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

				for (int i = 0; i < 3; i++) {
					SurfaceSmoothing::smooth(V, F, vertex_boundary_markers, part_of_original_stroke, new_mapped_indices, sharp_edge, dirty_boundary);
				}

				//Update the stroke positions after smoothing, in case their positions have changed (although they really shouldn't)
				initial_stroke->update_Positions(V);
				for (int i = 0; i < stroke_collection.size() - 1; i++) {
					stroke_collection[i].update_Positions(V);
				}

				viewervr.data.clear_all();
				viewervr.data.set_mesh_with_floor(V, F);

				draw_all_strokes();

				extrusion_base_already_drawn = false; //Reset
			}
			else { //mouse released after extrusion base drawn
				if (!extrusion_base->has_points_on_mesh) {
					hand_has_moved = false;
					viewervr.draw_while_computing = false;
					return true;
				}
				//extrusion_base->resample_all(); //This will shrink the drawn stroke. Might result in no face being contained inside the stroke
				dirty_boundary = true;
				extrusion_base->toLoop();
				extrusion_base_already_drawn = true;
				MeshExtrusion::extrude_prepare(*extrusion_base, base_surface_path); //Don't need to update all strokes here, since it didn't remove any vertices

				//base_model = extrusion_base->viewervr.corevr.model;
				//base_view = extrusion_base->viewervr.corevr.view;
				//base_proj = extrusion_base->viewervr.corevr.proj;
				base_model = extrusion_base->viewervr.corevr.get_model();
				base_view = extrusion_base->viewervr.corevr.get_view();
				base_proj = extrusion_base->viewervr.corevr.get_proj();
				base_viewport = extrusion_base->viewervr.corevr.viewport;

				draw_all_strokes();
				draw_extrusion_base(); //Need to draw the extrusion base separately, since it isn't added to the stroke_collection yet.

			}
		}
		else if (prev_tool_mode == TOGGLE || prev_tool_mode == FAIL) {
			viewervr.draw_while_computing = false;
		}


		prev_tool_mode = NONE;
		viewervr.draw_while_computing = false;
		cout << "draw while comp set to false" << endl;
		return true;
	}



	return true;
}


//TODO: make callback for this in viewer, like in exercise 5 of shapemod
bool callback_load_mesh(ViewerVR& viewervr, string filename) {
	/*	igl::readOFF(filename, V, F);
		viewervr.data.clear();
		viewervr.data.set_mesh(V, F);
		viewervr.data.compute_normals();
		viewervr.core.align_camera_center(viewervr.data.V);

		std::cout << filename.substr(filename.find_last_of("/") + 1) << endl;*/
	return true;
}

int main(int argc, char *argv[]) {
	//ViewerVR viewervr;
	//igl::read_triangle_mesh("../data/cube.off", V, F);
	viewervr.data.set_mesh_with_floor(V, F);

	viewervr.corevr.point_size = 55;
	//viewer.callback_load_mesh = callback_load_mesh;



	//Init stroke selector
	initial_stroke = new Stroke(V, F, viewervr, 0);
	if (argc == 2) {
		// Read mesh
	//	igl::readOFF(argv[1], V, F);
		//	callback_load_mesh(viewer, argv[1]);
	}
	else {
		// Read mesh
		//callback_load_mesh(viewer, "../data/cube.off");
	}

	//	callback_key_down(viewervr, '1', 0);

	CurveDeformation::smooth_deform_mode = true;
	viewervr.init();
	viewervr.callback_button_down = button_down;

	viewervr.launch();
}

