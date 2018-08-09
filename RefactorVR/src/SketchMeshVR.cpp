#ifdef _WIN32
#define NOMINMAX
#include <Windows.h>
#else
#include <unistd.h>
#include <curses.h>
#endif

#define STB_IMAGE_IMPLEMENTATION
#include "igl/opengl/glfw/imgui/stb_image.h"

#include <igl/readOFF.h>
#include <igl/edge_topology.h>
#include <igl/cat.h>
#include <igl/ray_mesh_intersect.h>
#include <iostream>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/per_corner_normals.h>

#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <igl/slice.h>
#include <SketchMeshVR.h>
#include <Stroke.h>
#include "SurfaceSmoothing.h"
#include "MeshCut.h"
#include "CurveDeformation.h"
#include "MeshExtrusion.h"
#include <igl/opengl/oculusVR.h>
#include "Mesh.h"
#include "Patch.h"
#include "LaplacianRemesh.h"

using namespace std;
using Viewer = igl::opengl::glfw::Viewer;
using OculusVR = igl::opengl::OculusVR;

OculusVR tmp;
Viewer viewer(tmp);
igl::opengl::glfw::imgui::ImGuiMenu menu;

// Vertex array, #V x3
Eigen::MatrixXd V(0, 3);
// Face array, #F x3
Eigen::MatrixXi F(0, 3);;
//Per vertex indicator of whether vertex position is fixed
Eigen::VectorXi vertex_is_fixed;
//Per edge indicator of what curve (or none if == 0) each edge belongs to
Eigen::VectorXi edge_boundary_markers;
//Per vertex indicator of whether vertex is on original stroke (outline of shape) (on OG stroke if ==1)
//Eigen::VectorXi part_of_original_stroke;
//Per edge indicator of whether the edge is sharp (if == 1 then sharp, otherwise smooth)
Eigen::VectorXi sharp_edge;
//Takes care of index mapping from before a cut/extrusion action to after (since some vertices are removed)
Eigen::VectorXi new_mapped_indices;

//Tracker needed for updating the strokes' closest vertex bindings after a stroke edge has been split (e.g. due to stroke add, cut or extrude)
Eigen::MatrixXi replacing_vertex_bindings(0, 4);

Mesh* base_mesh;

std::function<void()> current_tool_HUD;
std::function<void(void)> menu_HUD;
//General
enum ToolMode { DRAW, ADD, CUT, EXTRUDE, PULL, REMOVE, CHANGE, SMOOTH, NAVIGATE, NONE, DEFAULT, FAIL }; //NONE is used to indicate that a button was released, whereas DEFAULT indicates that one of the toggle buttons was pressed within its cooldown period, FAIL is used to indicate that something went wrong (e.g. user clicked too far away for PULL)

ToolMode tool_mode = DRAW, selected_tool_mode = DRAW, prev_tool_mode = NONE;
//Stroke* initial_stroke;
Stroke* added_stroke;
Stroke* extrusion_base;
vector<Stroke> stroke_collection;

//For smoothing
int initial_smooth_iter = 30;

//For selecting vertices
int handleID = -1;

//Variables for pulling a curve (and removing added control curves)
bool dirty_boundary = false;
int closest_stroke_ID, prev_closest_stroke_ID, prev_closest_stroke_idx = -1;

//Keeps track of the stroke IDs
int next_added_stroke_ID = 2; //Start at 2 because marker 1 belongs to the original boundary

//Variables for removing a control curve
bool stroke_was_removed = false;
int remove_stroke_clicked = 0;
Eigen::RowVector3d old_selected_stroke_color;

//Variables for cutting
bool cut_stroke_already_drawn = false;

//Variables for extrusion
bool extrusion_base_already_drawn = false;
SurfacePath base_surface_path;
Eigen::Matrix4f base_model, base_view, base_proj;
Eigen::Vector4f base_viewport;

//std::chrono::steady_clock::time_point _start_time, _end_time;
bool has_recentered = false;
bool draw_should_block = false;

Eigen::RowVector3d red(1, 0, 0);
Eigen::RowVector3d black(0, 0, 0);
Eigen::RowVector3d blue(0, 0, 1);

Eigen::RowVector3d laser_end_point;
bool prev_laser_show;

void draw_all_strokes() {
	viewer.data().set_edges(Eigen::MatrixXd(), Eigen::MatrixXi(), red); //Clear the non-original stroke edges
	Eigen::MatrixXd added_points;
	int points_to_hold_back;
	for (int i = 0; i < stroke_collection.size(); i++) {
		Eigen::VectorXi col_slice(3);
		col_slice.col(0) << 0, 1, 2;
		Eigen::MatrixXd edges_start, edges_end;
		Eigen::MatrixXi stroke_edges = stroke_collection[i].get_stroke_edges();
		igl::slice(V, stroke_edges.col(0).topRows(stroke_edges.rows() - !stroke_collection[i].is_loop), col_slice, edges_start);
		igl::slice(V, stroke_edges.col(1).topRows(stroke_edges.rows() - !stroke_collection[i].is_loop), col_slice, edges_end);
		//std::cout << "starts of: " << i << "    "<< stroke_edges.col(0).topRows(stroke_edges.rows() - !stroke_collection[i].is_loop).transpose() << std::endl;
		//std::cout << "ends of: " << i << "    " << stroke_edges.col(1).topRows(stroke_edges.rows() - !stroke_collection[i].is_loop).transpose() << std::endl;

		viewer.data().add_edges(edges_start, edges_end, stroke_collection[i].stroke_color);
	}
}

void draw_extrusion_base() {
	std::cout << "extru base" << std::endl;
	//The extrusion base does not get resampled, so use simple old method
	Eigen::MatrixXd added_points = extrusion_base->get3DPoints();
	std::cout << added_points << std::endl;
	int points_to_hold_back = 1 + !extrusion_base->is_loop;
	viewer.data().add_edges(added_points.topRows(added_points.rows() - points_to_hold_back), added_points.middleRows(1, added_points.rows() - points_to_hold_back), red);
}

void select_dragging_handle(Eigen::Vector3f& pos) {
	double closest_dist = INFINITY;
	double current_closest = closest_dist;
	closest_stroke_ID = -1;
	handleID = -1;
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

void set_laser_points(Eigen::Vector3f& pos) {
	Eigen::MatrixX3d LP(2, 3);
	Eigen::MatrixXd laser_color(2, 3);
	vector<igl::Hit> hits;

	if (igl::ray_mesh_intersect(pos, viewer.oculusVR.get_right_touch_direction(), V, F, hits)) { //Intersect the ray from the Touch controller with the mesh to get the 3D point
		laser_end_point = (V.row(F(hits[0].id, 0))*(1.0 - hits[0].u - hits[0].v) + V.row(F(hits[0].id, 1))*hits[0].u + V.row(F(hits[0].id, 2))*hits[0].v);
	}
	else { //First check for intersections with the mesh, then with the floor and finally just set the end point at a far distance_to_vert
		viewer.selected_data_index = 0;
		if (igl::ray_mesh_intersect(pos, viewer.oculusVR.get_right_touch_direction(), viewer.data().V, viewer.data().F, hits)) {
			laser_end_point = (viewer.data().V.row(viewer.data().F(hits[0].id, 0))*(1.0 - hits[0].u - hits[0].v) + viewer.data().V.row(viewer.data().F(hits[0].id, 1))*hits[0].u + viewer.data().V.row(viewer.data().F(hits[0].id, 2))*hits[0].v);
		}
		else {
			laser_end_point = (pos + 10 * viewer.oculusVR.get_right_touch_direction()).cast<double>();
		}
		viewer.selected_data_index = 1;
	}
	viewer.selected_data_index = 2;
	LP.row(0) = pos.cast<double>();
	LP.row(1) = laser_end_point;
	laser_color.setZero();
	viewer.data().set_laser_points(LP, laser_color);
	viewer.selected_data_index = 1;
}

void reset_before_draw() {
	cut_stroke_already_drawn = false;
	extrusion_base_already_drawn = false;
	(*base_mesh).V.resize(0, 3);
	(*base_mesh).F.resize(0, 3);
	(*base_mesh).vertex_is_fixed.resize(0);
	(*base_mesh).edge_boundary_markers.resize(0);
	(*base_mesh).new_mapped_indices.resize(0);
	(*base_mesh).sharp_edge.resize(0);
	(*base_mesh).mesh_to_patch_indices.resize(0);
	(*base_mesh).patches.clear(); //Request new patches
	(*base_mesh).face_patch_map.clear();
	viewer.data().clear();
	added_stroke->stroke_color = blue;
}

void sound_error_beep() {
#ifdef _WIN32
	Beep(500, 200);
#else
	beep();
#endif		
}

void button_down(OculusVR::ButtonCombo pressed, Eigen::Vector3f& pos) {
	set_laser_points(pos);
	if (pressed == OculusVR::ButtonCombo::TRIG && (selected_tool_mode == ADD || selected_tool_mode == REMOVE || selected_tool_mode == CUT || selected_tool_mode == EXTRUDE)) {
		if (stroke_collection.size() == 0) { //Don't go into these modes when there is no mesh yet
			prev_tool_mode = FAIL;
			return;
		}
	}
	else if (pressed == OculusVR::ButtonCombo::TRIG && selected_tool_mode == REMOVE) {
		if (stroke_collection.size() == 1) {
			return;
		}
	}
	else if (pressed == OculusVR::ButtonCombo::THUMB_MOVE) {
		V = viewer.data().V;
		for (int i = 0; i < stroke_collection.size(); i++) {
			stroke_collection[i].update_Positions(V, false);
		}
		draw_all_strokes();
		prev_tool_mode = NAVIGATE;
		return;
	}

	if (pressed == OculusVR::ButtonCombo::TRIG) {
		tool_mode = selected_tool_mode;
	}
	else if (pressed == OculusVR::ButtonCombo::NONE) {
		tool_mode = NONE;
	}
	else if (pressed == OculusVR::ButtonCombo::B) {
		tool_mode = SMOOTH;
	}


	if (tool_mode == DRAW) {
		if (draw_should_block) { //User has been too close to first sample point (closing the stroke too much), so we're in blocked state till the buttons are released again
			return;
		}
		if (prev_tool_mode == NONE) {
			if (has_recentered) {
				viewer.oculusVR.set_start_action_view(viewer.core.get_view());
				stroke_collection.clear();	
				added_stroke = new Stroke(V, F, viewer, 1);
				added_stroke->stroke_color = blue;
				next_added_stroke_ID = 2;
				added_stroke->addSegment(pos);
				prev_tool_mode = DRAW;
			}
			else {
				reset_before_draw();
				viewer.oculusVR.request_recenter();
				has_recentered = true;
			}
		}
		else if (prev_tool_mode == DRAW) {
			draw_should_block = added_stroke->addSegment(pos);
			return;
		}
	}
	else if (tool_mode == ADD) {
		if (prev_tool_mode == NONE || prev_tool_mode == FAIL) { //Adding a new control curve onto an existing mesh
			added_stroke = new Stroke(V, F, viewer, next_added_stroke_ID);
			added_stroke->stroke_color = blue;
			next_added_stroke_ID++;
			added_stroke->is_loop = false;
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

			if (igl::ray_mesh_intersect(pos, viewer.oculusVR.get_right_touch_direction(), V, F, hits)) { //Intersect the ray from the Touch controller with the mesh to get the 3D point
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
			for (int i = 1; i < stroke_collection.size(); i++) { //Skip stroke 0
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
			else if (closest_stroke_ID == prev_closest_stroke_ID) {
				remove_stroke_clicked++;
				draw_all_strokes();
			}
			else {
				remove_stroke_clicked = 1; //Start from 1
				if (prev_closest_stroke_idx != -1) {
					stroke_collection[prev_closest_stroke_idx].stroke_color = old_selected_stroke_color; //Reset the color of the previously selected stroke
				}
				//Redraw the original stroke and all added strokes, where the selected stroke is drawn in black.
				old_selected_stroke_color = stroke_collection[closest_stroke_idx].stroke_color;
				stroke_collection[closest_stroke_idx].stroke_color = black;
				prev_closest_stroke_idx = closest_stroke_idx;
				prev_closest_stroke_ID = closest_stroke_ID;
				draw_all_strokes();
			}

			if (remove_stroke_clicked == 2) { //Mechanism to force the user to click twice on the same stroke before removing it (safeguard)
				stroke_was_removed = true;
				stroke_collection[closest_stroke_idx].undo_stroke_add(edge_boundary_markers, sharp_edge, vertex_is_fixed);
	
				//We might have changed the patch structure (e.g. when removing a sharp cut stroke), so request new patches
				(*base_mesh).patches.clear();
				(*base_mesh).face_patch_map.clear();
				(*base_mesh).patches = Patch::init_patches(*base_mesh);

				stroke_collection.erase(stroke_collection.begin() + closest_stroke_idx);
				remove_stroke_clicked = 0; //Reset
				prev_closest_stroke_idx = -1;
			}

			prev_tool_mode = REMOVE;
		}
		if (prev_tool_mode == REMOVE) {
			return; //For REMOVE we only take action upon button press and release 
		}
	}
	else if (tool_mode == PULL) { //Dragging an existing curve
		//FAIL below is needed when initialy hand position is too far away from the curves. Switch to FAIL instead of NONE to ensure proper button-release handling
		if (prev_tool_mode == NONE || prev_tool_mode == FAIL) { //Also allow to go to pull after ADD because sometimes the buttons are hard to differentiate
			viewer.data().set_points(Eigen::MatrixXd(), black); //Will remove any previous unfinished strokes (e.g. unfinished cut)
			draw_all_strokes();
			select_dragging_handle(pos);

			if (handleID == -1) {//User clicked too far from any of the stroke vertices
				prev_tool_mode = FAIL;
				return;
			}
			CurveDeformation::startPullCurve(handleID, (*base_mesh).V, (*base_mesh).F);
			prev_tool_mode = PULL;
		}
		else if (prev_tool_mode == PULL) {
			CurveDeformation::pullCurve(pos.transpose().cast<double>(), (*base_mesh).V, (*base_mesh).edge_boundary_markers);
			for (int i = 0; i < (*base_mesh).patches.size(); i++) {
				(*base_mesh).patches[i]->update_patch_vertex_positions((*base_mesh).V);
			}

			for (int i = 0; i < 18; i++) {
				SurfaceSmoothing::smooth(*base_mesh, dirty_boundary);
			}

			//initial_stroke->update_Positions(V);
			for (int i = 0; i < stroke_collection.size(); i++) {
				stroke_collection[i].update_Positions(V, false);
			}

			viewer.data().set_mesh(V, F);
			Eigen::MatrixXd N_corners;
			igl::per_corner_normals(V, F, 50, N_corners);
			viewer.data().set_normals(N_corners);
			draw_all_strokes();
		}
	}
	else if (tool_mode == CUT) {
		if (prev_tool_mode == NONE || prev_tool_mode == FAIL) {
			prev_tool_mode = CUT;
			if (cut_stroke_already_drawn) { //clicked while cut stroke already drawn
				return;
			}

			draw_all_strokes(); //Will remove a possible old black (wrong) stroke
			added_stroke = new Stroke(V, F, viewer, next_added_stroke_ID); //For CUT, only increase the next added stroke ID upon button release (when we know it succeeded)
			added_stroke->stroke_color = red;
			viewer.oculusVR.set_start_action_view(viewer.core.get_view());
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
			if (extrusion_base_already_drawn) {
				added_stroke = new Stroke(V, F, viewer, next_added_stroke_ID);
				added_stroke->stroke_color = blue;
				next_added_stroke_ID++;
				added_stroke->addSegmentExtrusionSilhouette(pos);
			}
			else {
				extrusion_base = new Stroke(V, F, viewer, next_added_stroke_ID);
				extrusion_base->stroke_color = red;
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

			//if (initial_stroke->toLoop()) { //Returns false if the stroke only consists of 1 point (user just clicked)
			//	initial_stroke->generate3DMeshFromStroke(edge_boundary_markers, vertex_is_fixed, V, F);
			if (added_stroke->toLoop()) { //Returns false if the stroke only consists of 1 point (user just clicked)
				added_stroke->generate3DMeshFromStroke(edge_boundary_markers, vertex_is_fixed, V, F);

				if (!igl::is_edge_manifold(F)) { //Check if the drawn stroke results in an edge-manifold mesh, otherwise sound a beep and revert
					sound_error_beep();
				//	Eigen::MatrixXd drawn_points = initial_stroke->get3DPoints();
				//	initial_stroke->strokeReset();
					Eigen::MatrixXd drawn_points = added_stroke->get3DPoints();
				//	added_stroke->strokeReset();
					vertex_is_fixed.resize(0);
					edge_boundary_markers.resize(0);
					dirty_boundary = true;

					viewer.data().add_edges(drawn_points.block(0, 0, drawn_points.rows() - 1, 3), drawn_points.block(1, 0, drawn_points.rows() - 1, 3), black); //Display the stroke in black to show that it went wrong
					prev_tool_mode = NONE;
					std::cerr << "This stroke results in a mesh that is not edge manifold. Will not be able to process, try again (don't use self-intersecting strokes)." << std::endl;
					viewer.update_screen_while_computing = false;
					return;
				}

				Eigen::MatrixXi EV, FE, EF;
				igl::edge_topology(V, F, EV, FE, EF);
				sharp_edge.resize(EV.rows());
				sharp_edge.setZero(); //Set all edges to smooth after initial draw

				stroke_collection.push_back(*added_stroke);
				(*base_mesh).patches = Patch::init_patches(*base_mesh);

				dirty_boundary = true;
				for (int i = 0; i < initial_smooth_iter; i++) {
					SurfaceSmoothing::smooth(*base_mesh, dirty_boundary);
				}

				stroke_collection.back().update_Positions(V, true);
				//initial_stroke->update_Positions(V);
				viewer.data().clear();
				viewer.data().set_mesh(V, F);
				Eigen::MatrixXd N_corners;
				igl::per_corner_normals(V, F, 50, N_corners);
				viewer.data().set_normals(N_corners);

				//Overlay the drawn stroke
				//Eigen::MatrixXd added_points = initial_stroke->get3DPoints();
				//int points_to_hold_back = 1 + !initial_stroke->is_loop;
			//	Eigen::MatrixXd added_points = added_stroke->get3DPoints();
			//	int points_to_hold_back = 1 + !added_stroke->is_loop;
				//viewer.data().add_edges(added_points.topRows(added_points.rows() - points_to_hold_back), added_points.middleRows(1, added_points.rows() - points_to_hold_back), red);
				/*Eigen::VectorXi col_slice(3);
				col_slice.col(0) << 0, 1, 2;
				Eigen::MatrixXd edges_start, edges_end;
				Eigen::MatrixXi stroke_edges = added_stroke->get_stroke_edges();
				igl::slice(V, stroke_edges.col(0).topRows(stroke_edges.rows() - !added_stroke->is_loop), col_slice, edges_start);
				igl::slice(V, stroke_edges.col(1).topRows(stroke_edges.rows() - !added_stroke->is_loop), col_slice, edges_end);

				viewer.data().add_edges(edges_start, edges_end, red);*/

				draw_all_strokes();
			}
		}
		else if (prev_tool_mode == ADD) {
			dirty_boundary = true;
			if (!added_stroke->has_points_on_mesh) {
				sound_error_beep();
				std::cerr << "The stroke you are trying to add does not have any points that lie on the mesh surface. Will not be able to process, try again." << std::endl;
				viewer.update_screen_while_computing = false;
				return;
			}

			bool success;
			if (added_stroke->starts_on_mesh && added_stroke->ends_on_mesh) {
				success = LaplacianRemesh::remesh_open_path(*base_mesh, *added_stroke, replacing_vertex_bindings);
			}
			else if (!added_stroke->starts_on_mesh && !added_stroke->ends_on_mesh) { //Stroke like a cut stroke (starts and ends off mesh to wrap around)
			   //Need to remesh like it's a cut but without removing the inside faces (we get 2 loops of boundary vertices that both need to be stitched, so can't use remesh_open_path)
				added_stroke->prepend_first_point();
				added_stroke->append_final_point();
				added_stroke->toLoop();

				success = LaplacianRemesh::remesh_cutting_path(*base_mesh, *added_stroke, replacing_vertex_bindings);
			}
			else {
				std::cerr << "An added stroke either needs both the end- & startpoint to be outside of the mesh, or both on the mesh. Please try again. " << std::endl;
				success = false;
			}

			if (!success) {
				next_added_stroke_ID--; //Undo ID increment since stroke didn't actually get pushed back
				prev_tool_mode = NONE;
				sound_error_beep();
				viewer.data().set_points(Eigen::MatrixXd(), black); //Will remove the drawn stroke
				draw_all_strokes();
				Eigen::MatrixXd drawn_points = added_stroke->get3DPoints();
				int nr_edges = drawn_points.rows() - 1;
				viewer.data().add_edges(drawn_points.topRows(nr_edges), drawn_points.middleRows(1, nr_edges), black); //Display the stroke in black to show that it went wrong
				viewer.update_screen_while_computing = false;
				return;
			}

			stroke_collection.push_back(*added_stroke);
			(*base_mesh).patches.clear();
			(*base_mesh).face_patch_map.clear();
			(*base_mesh).patches = Patch::init_patches(*base_mesh);

			//initial_stroke->update_vert_bindings(new_mapped_indices, edge_boundary_markers, sharp_edge, vertex_is_fixed, replacing_vertex_bindings); //Don't test if the initial one dies

			int nr_removed = 0, original_collection_size = stroke_collection.size();
			for (int i = 0; i < original_collection_size - 1; i++) { //Don't update the added stroke, as this is done inside the remeshing already
				if (!stroke_collection[i - nr_removed].update_vert_bindings(new_mapped_indices, edge_boundary_markers, sharp_edge, vertex_is_fixed, replacing_vertex_bindings)) { //Stroke dies, don't need to do stroke.undo_stroke_add, cause all its vertices also cease to exist
					stroke_collection.erase(stroke_collection.begin() + i - nr_removed);
					nr_removed++;
					dirty_boundary = true;
				}
			}

			//Update the stroke positions, will update with resampled stroke points
		//	initial_stroke->update_Positions(V);
			for (int i = 0; i < stroke_collection.size(); i++) { //TODO: I removed -1 here, check if still works
				stroke_collection[i].update_Positions(V, true);
			}

			viewer.data().clear();
			viewer.data().set_mesh(V, F);
			Eigen::MatrixXd N_corners;
			igl::per_corner_normals(V, F, 50, N_corners);
			viewer.data().set_normals(N_corners);

			draw_all_strokes();
		}
		else if (prev_tool_mode == REMOVE && stroke_was_removed) { //Only redraw if we actually removed a stroke (otherwise we draw unnecessary)
			dirty_boundary = true;
			for (int i = 0; i < 10; i++) {
				SurfaceSmoothing::smooth(*base_mesh, dirty_boundary);
			}

			for (int i = 0; i < stroke_collection.size(); i++) {
				stroke_collection[i].update_Positions(V, false);
			}

			//TODO: calling clear will erase the mesh/trackball rotation, so only do clear() when absolutely necessary. Check if below works like this
			//viewer.data().clear();
			viewer.data().set_mesh(V, F);
			Eigen::MatrixXd N_corners;
			igl::per_corner_normals(V, F, 50, N_corners);
			viewer.data().set_normals(N_corners);

			stroke_was_removed = false;
			draw_all_strokes();
		}
		else if (prev_tool_mode == PULL && handleID != -1) {
			for (int i = 0; i < 20; i++) {
				SurfaceSmoothing::smooth(*base_mesh, dirty_boundary);
			}

		//	initial_stroke->update_Positions(V);
			for (int i = 0; i < stroke_collection.size(); i++) {
				stroke_collection[i].update_Positions(V, false);
			}

			//viewer.data().clear();
			viewer.data().set_mesh(V, F);
			Eigen::MatrixXd N_corners;
			igl::per_corner_normals(V, F, 50, N_corners);
			viewer.data().set_normals(N_corners);
			draw_all_strokes();
		}
		else if (prev_tool_mode == CUT) {
			if (!added_stroke->has_points_on_mesh) {
				prev_tool_mode = NONE;
				viewer.update_screen_while_computing = false;
				return;
			}
			if (cut_stroke_already_drawn) { //User had already drawn the cut stroke and has now clicked/drawn the final stroke for removing the part. User should click on the part that should be removed
				dirty_boundary = true;
				vector<igl::Hit> hits;
				if (!igl::ray_mesh_intersect(pos, viewer.oculusVR.get_right_touch_direction(), V, F, hits)) { //Intersect the ray from the Touch controller with the mesh to get the 3D point
					viewer.update_screen_while_computing = false;
					return;
				}
				int clicked_face = hits[0].id;

				added_stroke->prepend_first_point();
				added_stroke->append_final_point();
				added_stroke->toLoop();

				bool cut_success = MeshCut::cut(*base_mesh, *added_stroke, clicked_face, replacing_vertex_bindings);
				if (!cut_success) { //Catches the following cases: when the cut removes all mesh vertices/faces, when the first/last cut point aren't correct (face -1 in SurfacePath) or when sorting takes "infinite" time
					cut_stroke_already_drawn = false;
					prev_tool_mode = NONE;
					sound_error_beep();
					viewer.data().set_points(Eigen::MatrixXd(), black); //Will remove the pink cut stroke
					draw_all_strokes();
					Eigen::MatrixXd drawn_points = added_stroke->get3DPoints();
					int nr_edges = drawn_points.rows() - 1;
					viewer.data().add_edges(drawn_points.topRows(nr_edges), drawn_points.middleRows(1, nr_edges), black); //Display the stroke in black to show that it went wrong
					viewer.update_screen_while_computing = false;
					return;
				}
				(*base_mesh).patches.clear();
				(*base_mesh).face_patch_map.clear();
				(*base_mesh).patches = Patch::init_patches(*base_mesh);

				stroke_collection.push_back(*added_stroke);
				next_added_stroke_ID++; //Stroke is successfully pushed back

			//	initial_stroke->update_vert_bindings(new_mapped_indices, edge_boundary_markers, sharp_edge, vertex_is_fixed, replacing_vertex_bindings); //Don't test if the initial one dies

				int nr_removed = 0, original_collection_size = stroke_collection.size();
				for (int i = 0; i < original_collection_size; i++) {
					if (!stroke_collection[i - nr_removed].update_vert_bindings(new_mapped_indices, edge_boundary_markers, sharp_edge, vertex_is_fixed, replacing_vertex_bindings)) { //Stroke dies, don't need to do stroke.undo_stroke_add, cause all its vertices also cease to exist
						stroke_collection.erase(stroke_collection.begin() + i - nr_removed);
						nr_removed++;
						dirty_boundary = true;
					}
				}

				for (int i = 0; i < 10; i++) {
					std::cout << "iter: " << i << std::endl;
					SurfaceSmoothing::smooth(*base_mesh, dirty_boundary);
				}

				//Update the stroke positions after smoothing, will also add resampled points
				//initial_stroke->update_Positions(V);
				for (int i = 0; i < stroke_collection.size(); i++) {
					stroke_collection[i].update_Positions(V, true);
				}

				viewer.data().clear();
				viewer.data().set_mesh(V, F);
				Eigen::MatrixXd N_corners;
				igl::per_corner_normals(V, F, 50, N_corners);
				viewer.data().set_normals(N_corners);

				cut_stroke_already_drawn = false;
				draw_all_strokes();
			}
			else { //We're finished drawing the cut stroke, prepare for when user draws the final stroke to remove the part
				if (added_stroke->has_self_intersection()) {
					std::cerr << "Cut stroke contains a loop which is illegal. Please try again. " << std::endl;
					prev_tool_mode = NONE;
					sound_error_beep();
					viewer.data().set_points(Eigen::MatrixXd(), black); //Will remove the pink cut stroke
					draw_all_strokes();
					Eigen::MatrixXd drawn_points = added_stroke->get3DPoints();
					int nr_edges = drawn_points.rows() - 1;
					viewer.data().add_edges(drawn_points.topRows(nr_edges), drawn_points.middleRows(1, nr_edges), black); //Display the stroke in black to show that it went wrong
					viewer.update_screen_while_computing = false;
					return;
				}
				cut_stroke_already_drawn = true;
			}
		}
		else if (prev_tool_mode == EXTRUDE) {
			if (extrusion_base_already_drawn) { //User has drawn the silhouette stroke for extrusion
				dirty_boundary = true;				
				added_stroke->toLoop();
				added_stroke->is_loop = false; //Set to false manually, because we don't want curveDeformation to consider it as a loop (but we do need looped 3DPoints)

				bool success_extrude = false;
				if (added_stroke->get3DPoints().rows() > 6) {
					success_extrude = MeshExtrusion::extrude_main(*base_mesh, base_surface_path, *added_stroke, *extrusion_base, base_model, base_view, base_proj, base_viewport, replacing_vertex_bindings);
				}

				if (!success_extrude) { //Catches the case that the extrusion base removes all faces/vertices or that the extrusion silhouette has too little samples/is just a click
					viewer.update_screen_while_computing = false;
					prev_tool_mode = NONE;
					next_added_stroke_ID -= 2;
					extrusion_base_already_drawn = false;
					sound_error_beep();
					draw_all_strokes(); //Will remove the drawn base & silhouette strokes
					Eigen::MatrixXd drawn_points = extrusion_base->get3DPoints();
					int nr_edges = drawn_points.rows() - 1;
					viewer.data().add_edges(drawn_points.topRows(nr_edges), drawn_points.middleRows(1, nr_edges), black); //Display the stroke in black to show that it went wrong
					drawn_points = added_stroke->get3DPoints();
					nr_edges = drawn_points.rows() - 1;
					viewer.data().add_edges(drawn_points.topRows(nr_edges - 1), drawn_points.middleRows(1, nr_edges - 1), black); //Display the stroke in black to show that it went wrong
					viewer.selected_data_index = 2;
					viewer.data().show_laser = true;
					viewer.selected_data_index = 1;
					return;
				}
				(*base_mesh).patches.clear();
				(*base_mesh).face_patch_map.clear();
				(*base_mesh).patches = Patch::init_patches(*base_mesh);

				stroke_collection.push_back(*extrusion_base);
				stroke_collection.push_back(*added_stroke);

				//initial_stroke->update_vert_bindings(new_mapped_indices, edge_boundary_markers, sharp_edge, vertex_is_fixed, replacing_vertex_bindings);

				int nr_removed = 0, original_collection_size = stroke_collection.size();
				for (int i = 0; i < original_collection_size - 1; i++) { //Skip the newly added stroke since it is already updated inside extrude_main()
					if (!stroke_collection[i - nr_removed].update_vert_bindings(new_mapped_indices, edge_boundary_markers, sharp_edge, vertex_is_fixed, replacing_vertex_bindings)) {
						//Stroke dies, don't need to do stroke.undo_stroke_add, cause all its vertices also cease to exist (or in case of non-loop strokes that have a middle portion removed, the undo_stroke_add is done inside of the update_vert_bindings)
						stroke_collection.erase(stroke_collection.begin() + i - nr_removed);
						nr_removed++;
						dirty_boundary = true;
					}
				}

				for (int i = 0; i < 10; i++) {
						SurfaceSmoothing::smooth(*base_mesh, dirty_boundary);
				}

				//Update the stroke positions after smoothing, also updates resampled stroke points
				//initial_stroke->update_Positions(V);
				for (int i = 0; i < stroke_collection.size(); i++) { //TODO: I removed -1 here, check that it still works
					stroke_collection[i].update_Positions(V, true);
				}

				viewer.data().clear();
				viewer.data().set_mesh(V, F);
				Eigen::MatrixXd N_corners;
				igl::per_corner_normals(V, F, 50, N_corners);
				viewer.data().set_normals(N_corners);

				draw_all_strokes();
				viewer.selected_data_index = 2;
				viewer.data().show_laser = true;
				viewer.selected_data_index = 1;
				extrusion_base_already_drawn = false;
			}
			else { //mouse released after extrusion base drawn
				if (!extrusion_base->has_points_on_mesh) {
					viewer.update_screen_while_computing = false;
					prev_tool_mode = NONE;
					return;
				}
				else if (extrusion_base->has_been_outside_mesh || extrusion_base->has_self_intersection()) {
					if (extrusion_base->has_self_intersection()) {
						std::cerr << "Extrusion base contains a loop which is illegal. Please try again. " << std::endl;
					}
					else {
						std::cerr << "Extrusion base has parts which are drawn outside of the mesh surface, which is not allowed. Try again. " << std::endl;
					}
					prev_tool_mode = NONE;
					sound_error_beep();
					draw_all_strokes(); //Will remove the drawn base & silhouette strokes
					Eigen::MatrixXd drawn_points = extrusion_base->get3DPoints();
					int nr_edges = drawn_points.rows() - 1;
					viewer.data().add_edges(drawn_points.topRows(nr_edges), drawn_points.middleRows(1, nr_edges), black); //Display the stroke in black to show that it went wrong
					viewer.update_screen_while_computing = false;
					return;
				}

				dirty_boundary = true;
				extrusion_base->toLoop();

				bool succes_extrude_prepare = MeshExtrusion::extrude_prepare(*extrusion_base, base_surface_path); //Don't need to update all strokes here, since it didn't remove any vertices
				if (!succes_extrude_prepare) { //Catches the case that face == -1 in SurfacePath
					sound_error_beep();
					next_added_stroke_ID--;
					prev_tool_mode = NONE;
					draw_all_strokes(); //Removes the drawn base stroke
					Eigen::MatrixXd drawn_points = extrusion_base->get3DPoints();
					int nr_edges = drawn_points.rows();
					viewer.data().add_edges(drawn_points.topRows(nr_edges), drawn_points.middleRows(1, nr_edges), black); //Display the stroke in black to show that it went wrong
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
				viewer.selected_data_index = 2;
				viewer.data().show_laser = false;
				viewer.selected_data_index = 1;
			}
		}
		else if (prev_tool_mode == FAIL) {
			viewer.update_screen_while_computing = false;
		}
		else if (prev_tool_mode == SMOOTH) {
			if (V.rows() == 0) {
				prev_tool_mode = NONE;
				viewer.update_screen_while_computing = false;
				return;
			}
			SurfaceSmoothing::smooth(*base_mesh, dirty_boundary);
			for (int i = 0; i < stroke_collection.size(); i++) {
				stroke_collection[i].update_Positions(V, false);
			}
			//	viewer.data().clear();
			viewer.data().set_mesh(V, F);
			Eigen::MatrixXd N_corners;
			igl::per_corner_normals(V, F, 50, N_corners);
			viewer.data().set_normals(N_corners);

			draw_all_strokes();
		}
		else if (prev_tool_mode == NAVIGATE) {
			V = viewer.data().V;
			//initial_stroke->update_Positions(V);
			for (int i = 0; i < stroke_collection.size(); i++) {
				stroke_collection[i].update_Positions(V, true);
			}	
			draw_all_strokes();

			if (extrusion_base_already_drawn) {
				//First need to update the extrusion base manually, since it doesn't have any closest_vertex_indices yet. Retrieve vertex positions
				vector<PathElement> base_path = base_surface_path.get_path();
				Eigen::MatrixX3d new_3DPoints(base_path.size(), 3);
				for (int i = 0; i < base_path.size(); i++) {
					new_3DPoints.row(i) = base_path[i].get_vertex().transpose();
				}
				extrusion_base->set3DPoints(new_3DPoints);
				draw_extrusion_base();
			}
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
	Eigen::MatrixXd N_corners;
	igl::per_corner_normals(V, F, 50, N_corners);
	viewer.data().set_normals(N_corners);

	std::cout << filename.substr(filename.find_last_of("/") + 1) << endl;
	return true;
}

void menu_opened() {
	menu.set_active();
	menu.callback_draw_viewer_window = menu_HUD;
	viewer.selected_data_index = 2;
	prev_laser_show = viewer.data().show_laser;
	viewer.data().show_laser = true;
	viewer.selected_data_index = 1;
	viewer.oculusVR.menu_active = true;
}

void menu_closed() {
	menu.set_inactive();
	menu.callback_draw_viewer_window = current_tool_HUD;
	viewer.selected_data_index = 2;
	viewer.data().show_laser = prev_laser_show;
	viewer.selected_data_index = 1;
	viewer.oculusVR.menu_active = false;
}

void reset_trackers() {
	if (extrusion_base_already_drawn) { //If we're switching from an unfinished extrusion
		next_added_stroke_ID--;
	}
	cut_stroke_already_drawn = false;
	extrusion_base_already_drawn = false;
	if (prev_closest_stroke_idx != -1) {
		stroke_collection[prev_closest_stroke_idx].stroke_color = old_selected_stroke_color; //Reset the color of the previously selected stroke
		prev_closest_stroke_ID = -1;
	}
	remove_stroke_clicked = 0;

	viewer.data().set_points(Eigen::MatrixXd(), black); //Will remove a possible pink cut stroke
	draw_all_strokes();
}

int main(int argc, char *argv[]) {
	//Init stroke selector
	added_stroke = new Stroke(V, F, viewer, 1);
	base_mesh = new Mesh(V, F, edge_boundary_markers, vertex_is_fixed, new_mapped_indices, sharp_edge, 0);

	Eigen::MatrixXd V_floor(4, 3);
	V_floor.row(0) << -10, 0, -10;
	V_floor.row(1) << 10, 0, -10;
	V_floor.row(2) << 10, 0, 10;
	V_floor.row(3) << -10, 0, 10;

	Eigen::MatrixXi F_floor(2, 3);
	F_floor.row(0) << 0, 3, 1;
	F_floor.row(1) << 3, 2, 1;

	char cur_dir[256];
	GetCurrentDirectoryA(256, cur_dir);

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
		V_uv = V_uv.array() * 2;

		Eigen::Matrix<unsigned char, Dynamic, Dynamic> texR, texG, texB;
		std::string texture_file = std::string(cur_dir) + "\\..\\data\\free\\floor.png";
		int width, height, n;
		unsigned char *data = stbi_load(texture_file.c_str(), &width, &height, &n, 4);

		if (!data) {
			std::cerr << "Could not load floor texture." << std::endl;
		}
		texR.resize(height, width);
		texG.resize(height, width);
		texB.resize(height, width);

		for (unsigned j = 0; j < height; ++j) {
			for (unsigned i = 0; i < width; ++i) {
				// used to flip with libPNG, but I'm not sure if
				// simply j*width + i wouldn't be better
				// stb_image uses horizontal scanline an starts top-left corner
				texR(i, j) = data[4 * ((width - 1 - i) + width * (height - 1 - j))];
				texG(i, j) = data[4 * ((width - 1 - i) + width * (height - 1 - j)) + 1];
				texB(i, j) = data[4 * ((width - 1 - i) + width * (height - 1 - j)) + 2];
			}
		}
		stbi_image_free(data);

		viewer.data().set_uv(V_uv);
		viewer.data().show_texture = true;
		viewer.data().set_texture(texR, texG, texB);
	}

	viewer.append_mesh();
	viewer.data().set_mesh(V, F);
	Eigen::MatrixXd N_corners;
	igl::per_corner_normals(V, F, 50, N_corners);
	viewer.data().set_normals(N_corners);
	viewer.append_mesh(); //For laser ray/point
	viewer.data().show_laser = false;
	viewer.selected_data_index = 1;
	viewer.plugins.push_back(&menu);

	//CurveDeformation::smooth_deform_mode = true;
	viewer.init_oculus();

	GLuint img_texture = 0, img_texture1 = 0, img_texture2 = 0, img_texture3 = 0, img_texture4 = 0, img_texture5 = 0, img_texture6 = 0;
	int img_width, img_height, nrChannels;
	std::string filename = std::string(cur_dir) + "\\..\\data\\free\\draw.png";
	unsigned char *img_data = stbi_load(filename.c_str(), &img_width, &img_height, &nrChannels, 4);
	if (!img_data) {
		std::cerr << "Could not load image 1." << std::endl;
	}

	glGenTextures(1, &img_texture);
	GLenum err = glGetError();
	glBindTexture(GL_TEXTURE_2D, img_texture);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, img_width, img_height, 0, GL_RGBA, GL_UNSIGNED_BYTE, img_data);
	void* im_texID_draw = (void *)(intptr_t)img_texture;

	filename = std::string(cur_dir) + "\\..\\data\\free\\pull.png";
	img_data = stbi_load(filename.c_str(), &img_width, &img_height, &nrChannels, 4);
	if (!img_data) {
		std::cerr << "Could not load image 2." << std::endl;
	}

	glGenTextures(1, &img_texture2);
	glBindTexture(GL_TEXTURE_2D, img_texture2);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, img_width, img_height, 0, GL_RGBA, GL_UNSIGNED_BYTE, img_data);
	void* im_texID_pull = (void *)(intptr_t)img_texture2;

	filename = std::string(cur_dir) + "\\..\\data\\free\\plus.png";
	img_data = stbi_load(filename.c_str(), &img_width, &img_height, &nrChannels, 4);
	if (!img_data) {
		std::cerr << "Could not load image 3." << std::endl;
	}

	glGenTextures(1, &img_texture3);
	glBindTexture(GL_TEXTURE_2D, img_texture3);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, img_width, img_height, 0, GL_RGBA, GL_UNSIGNED_BYTE, img_data);
	void* im_texID_add = (void *)(intptr_t)img_texture3;

	filename = std::string(cur_dir) + "\\..\\data\\free\\scissor.png";
	img_data = stbi_load(filename.c_str(), &img_width, &img_height, &nrChannels, 4);
	if (!img_data) {
		std::cerr << "Could not load image 4." << std::endl;
	}

	glGenTextures(1, &img_texture4);
	glBindTexture(GL_TEXTURE_2D, img_texture4);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, img_width, img_height, 0, GL_RGBA, GL_UNSIGNED_BYTE, img_data);
	void* im_texID_cut = (void *)(intptr_t)img_texture4;

	filename = std::string(cur_dir) + "\\..\\data\\free\\bump.png";
	img_data = stbi_load(filename.c_str(), &img_width, &img_height, &nrChannels, 4);
	if (!img_data) {
		std::cerr << "Could not load image 5." << std::endl;
	}

	glGenTextures(1, &img_texture5);
	glBindTexture(GL_TEXTURE_2D, img_texture5);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, img_width, img_height, 0, GL_RGBA, GL_UNSIGNED_BYTE, img_data);
	void* im_texID_extrude = (void *)(intptr_t)img_texture5;

	filename = std::string(cur_dir) + "\\..\\data\\free\\minus.png";
	img_data = stbi_load(filename.c_str(), &img_width, &img_height, &nrChannels, 4);
	if (!img_data) {
		std::cerr << "Could not load image 6." << std::endl;
	}

	glGenTextures(1, &img_texture6);
	glBindTexture(GL_TEXTURE_2D, img_texture6);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, img_width, img_height, 0, GL_RGBA, GL_UNSIGNED_BYTE, img_data);
	void* im_texID_remove = (void *)(intptr_t)img_texture6;

	void* im_texID_cur = im_texID_draw;

	current_tool_HUD = [&]() {
		ImGui::SetNextWindowSize(ImVec2(200.0f, 200.0f), ImGuiSetCond_FirstUseEver);
		ImGui::SetNextWindowPos(ImVec2(0.0f, 0.0f), ImGuiSetCond_FirstUseEver);
		ImGuiStyle& style = ImGui::GetStyle();
		style.WindowPadding = ImVec2(0, 0);
		style.WindowRounding = 0.0f;
		style.WindowTitleAlign = ImVec2(0.5f, 0.5f);
		style.ItemInnerSpacing = ImVec2(0, 0);
		style.Colors[ImGuiCol_WindowBg] = ImVec4(0.00f, 0.00f, 0.00f, 1.0f);
		ImGui::Begin("Current Tool", 0, ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoCollapse | ImGuiWindowFlags_NoScrollbar);
		ImGui::SameLine(0.5*(200 - 178)); //Take distance_to_vert around image and offset by half of it
		ImGui::Image(im_texID_cur, ImVec2(178.0f, 178.0f), ImVec2(0, 0), ImVec2(1, 1), ImColor(255, 255, 255, 255), ImColor(0, 0, 0, 0));
		ImGui::End();
	};

	menu_HUD = [=, &im_texID_cur]() {
		const ImVec2 texsize = ImGui_ImplGlfwGL3_GetTextureSize();
		ImColor icon_background_color(0, 0, 0, 255);
		ImVec2 uv_start(0, 0);
		ImVec2 uv_end(1, 1);
		ImVec2 icon_size(320, 320);
		//std::cout << "texsize" << texsize.x << "  " << texsize.y << std::endl;
		//ImGui::SetNextWindowSize(texsize);
		ImGuiStyle& style = ImGui::GetStyle();
		style.WindowPadding = ImVec2(0, 0);
		style.WindowRounding = 0.0f;
		style.WindowTitleAlign = ImVec2(0.5f, 0.5f);
		style.ItemInnerSpacing = ImVec2(0, 0);
		style.Colors[ImGuiCol_WindowBg] = ImVec4(0.00f, 0.00f, 0.00f, 1.0f);

		ImGui::SetNextWindowSize(ImVec2(1024.0f, 700.0f), ImGuiSetCond_FirstUseEver);
		ImGui::SetNextWindowPos(ImVec2(0.0f, 0.0f), ImGuiSetCond_FirstUseEver);

		ImGui::Begin("Selection Menu", 0, ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoCollapse | ImGuiWindowFlags_NoScrollbar);
		int frame_padding = 5;
		ImGuiIO& io = ImGui::GetIO();
		ImGui::SetWindowFontScale(2.f);

		ImGui::PushID(0);
		if (ImGui::ImageButton(im_texID_draw, icon_size, uv_start, uv_end, frame_padding, icon_background_color)) {
			reset_trackers();
			selected_tool_mode = DRAW;
			im_texID_cur = im_texID_draw;
			viewer.selected_data_index = 2;
			viewer.data().show_laser = false;
			prev_laser_show = viewer.data().show_laser;
			viewer.selected_data_index = 1;
			menu_closed();
		}
		ImGui::PopID();
		ImGui::SameLine();
		ImGui::PushID(1);
		//ImGui::SetCursorPos(ImVec2(0, 345));
		if (ImGui::ImageButton(im_texID_pull, icon_size, uv_start, uv_end, frame_padding, icon_background_color)) {
			reset_trackers();
			selected_tool_mode = PULL;
			im_texID_cur = im_texID_pull;
			viewer.selected_data_index = 2;
			viewer.data().show_laser = false;
			prev_laser_show = viewer.data().show_laser;
			viewer.selected_data_index = 1;
			menu_closed();
		}
		ImGui::PopID();
		ImGui::SameLine();
		ImGui::PushID(2);
		//ImGui::SetCursorPos(ImVec2(379, 379));
		if (ImGui::ImageButton(im_texID_add, icon_size, uv_start, uv_end, frame_padding, icon_background_color)) {
			reset_trackers();
			selected_tool_mode = ADD;
			im_texID_cur = im_texID_add;
			viewer.selected_data_index = 2;
			viewer.data().show_laser = true;
			prev_laser_show = viewer.data().show_laser;
			viewer.selected_data_index = 1;
			menu_closed();
		}
		ImGui::PopID();

		ImGui::PushID(3);
		if (ImGui::ImageButton(im_texID_cut, icon_size, uv_start, uv_end, frame_padding, icon_background_color)) {
			reset_trackers();
			selected_tool_mode = CUT;
			im_texID_cur = im_texID_cut;
			viewer.selected_data_index = 2;
			viewer.data().show_laser = true;
			prev_laser_show = viewer.data().show_laser;
			viewer.selected_data_index = 1;
			menu_closed();
		}
		ImGui::PopID();
		ImGui::SameLine();
		ImGui::PushID(4);
		//ImGui::SetCursorPos(ImVec2(0, 345));
		if (ImGui::ImageButton(im_texID_extrude, icon_size, uv_start, uv_end, frame_padding, icon_background_color)) {
			reset_trackers();
			selected_tool_mode = EXTRUDE;
			im_texID_cur = im_texID_extrude;
			viewer.selected_data_index = 2;
			viewer.data().show_laser = true;
			prev_laser_show = viewer.data().show_laser;
			viewer.selected_data_index = 1;
			menu_closed();
		}
		ImGui::PopID();
		ImGui::SameLine();
		ImGui::PushID(5);
		//ImGui::SetCursorPos(ImVec2(379, 379));
		if (ImGui::ImageButton(im_texID_remove, icon_size, uv_start, uv_end, frame_padding, icon_background_color)) {
			reset_trackers();
			selected_tool_mode = REMOVE;
			im_texID_cur = im_texID_remove;
			viewer.selected_data_index = 2;
			viewer.data().show_laser = true;
			prev_laser_show = viewer.data().show_laser;
			viewer.selected_data_index = 1;
			menu_closed();
		}
		ImGui::PopID();
		ImGui::End();
	};

	menu.callback_draw_viewer_window = current_tool_HUD;
	viewer.oculusVR.callback_button_down = button_down;
	viewer.oculusVR.callback_menu_opened = menu_opened;
	viewer.oculusVR.callback_menu_closed = menu_closed;
	viewer.data().point_size = 15;
	viewer.data().show_lines = true; //TODO change
	viewer.launch_oculus();
}

