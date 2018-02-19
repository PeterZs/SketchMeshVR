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
#include <iostream>
#include <unordered_map>
#include <igl/viewer/VR_Viewer.h>
#include <SketchMeshVR.h>
#include <Stroke.h>
#include "SurfaceSmoothing.h"


using namespace std;
using ViewerVR = igl::viewer::VR_Viewer;

// Vertex array, #V x3
Eigen::MatrixXd V;
// Face array, #F x3
Eigen::MatrixXi F;
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
enum ToolMode { DRAW, ADD, CUT, EXTRUDE, PULL, REMOVE, CHANGE, SMOOTH, NAVIGATE, NONE };

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

//Variables for adding control curves
bool last_add_on_mesh = false;
unordered_map<int, int> backside_vertex_map;

//Variables for removing a control curve
bool stroke_was_removed = false;
int remove_stroke_clicked = 0;

//Variables for cutting
bool cut_stroke_already_drawn = false;

void draw_all_strokes(ViewerVR& viewervr) {
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

bool button_down(ViewerVR::ButtonCombo pressed, Eigen::Vector3f& pos, igl::viewer::VR_Viewer& viewervr) {
    ToolMode pressed_type;
    if(pressed == ViewerVR::ButtonCombo::GRIPTRIG){
        pressed_type = DRAW;
    }else if(pressed == ViewerVR::ButtonCombo::GRIP){
		cout << "cutting" << endl;

		pressed_type = CUT;
    }else if(pressed == ViewerVR::ButtonCombo::TRIG){
		cout << "pulling" << endl;

		pressed_type = PULL;
    }else if(pressed == ViewerVR::ButtonCombo::A){
        
    }else if(pressed == ViewerVR::ButtonCombo::B){
        
    }else if(pressed == ViewerVR::ButtonCombo::X){
        
    }else if(pressed == ViewerVR::ButtonCombo::Y){
        
    }else if(pressed == ViewerVR::ButtonCombo::NONE){
        pressed_type = NONE;
    }

	if (!((pos - prev_pos).isZero())) {
		hand_has_moved = true;
	}
	prev_pos = pos;
    
	if (pressed_type == PULL || pressed_type == ADD || pressed_type == CUT || pressed_type == EXTRUDE) {
		if (initial_stroke->empty2D()) { //Don't go into these modes when there is no mesh yet
			return true;
		}
	}
	if (pressed_type == REMOVE) {
		if (stroke_collection.size() == 0) {
			return true;
		}
		remove_stroke_clicked = 0; //Reset because we might be left with a single click from the last round
	}
	else if (pressed_type == CUT) {
		cut_stroke_already_drawn = false; //Reset because we might have stopped before finishing the cut last time
	}
	tool_mode = pressed_type;

	if (tool_mode == DRAW) { //Creating the first curve/mesh
		if (prev_tool_mode == NONE) {
			viewervr.data.clear_without_floor();
			stroke_collection.clear();
			next_added_stroke_ID = 2;
			initial_stroke->strokeReset();
			initial_stroke->strokeAddSegment(pos);
			prev_tool_mode = DRAW;
			skip_standardcallback = true;
		}
		else if (prev_tool_mode == DRAW) {
			cout << "adding to stroke" << endl;
			//We had already started drawing, continue
			initial_stroke->strokeAddSegment(pos);
			return true;
		}
    }else if(tool_mode == ADD){
        if(prev_tool_mode == NONE){ //Adding a new control curve onto an existing mesh
            added_stroke = new Stroke(V, F, viewervr, next_added_stroke_ID);
            next_added_stroke_ID++;
            added_stroke->strokeAddSegmentAdd(pos); //If the user starts outside of the mesh, consider the movement as navigation
         skip_standardcallback = true;
        }else if(prev_tool_mode == ADD){
            last_add_on_mesh = added_stroke->strokeAddSegmentAdd(pos);
            return true;
        }
        
    }
	else if (tool_mode == NONE) {	//Have to finish up as if we're calling mouse_up()
		if (prev_tool_mode == NONE) {
			return true;
		}

		else if (prev_tool_mode == DRAW) {
			cout << "rounding up " << endl;
			initial_stroke->strokeAddSegment(pos);
			if (initial_stroke->toLoop()) {//Returns false if the stroke only consists of 1 point (user just clicked)
										   //Give some time to show the stroke
#ifdef _WIN32
				Sleep(200);
#else
				usleep(200000);
#endif
				cout << "made looop " << endl;
				backside_vertex_map = initial_stroke->generate3DMeshFromStroke(vertex_boundary_markers, part_of_original_stroke);
				F = viewervr.data.F.block(0,0,viewervr.data.F.rows()-2, viewervr.data.F.cols()); //Don't consider the last 2 faces because they belong to the floor
				V = viewervr.data.V.block(0, 0, viewervr.data.V.rows() - 4, viewervr.data.V.cols()); //Don't consider the last 4 vertices because they belong to the floor

				//TODO: enable this once mesh is a manifold again
			/*	Eigen::MatrixXi EV, FE, EF;
				igl::edge_topology(V, F, EV, FE, EF);
				sharp_edge.resize(EV.rows());
				sharp_edge.setZero(); //Set all edges to smooth after initial draw

				dirty_boundary = true;

				for (int i = 0; i < initial_smooth_iter; i++) {
					SurfaceSmoothing::smooth(V, F, vertex_boundary_markers, part_of_original_stroke, new_mapped_indices, sharp_edge, dirty_boundary);
				}
				*/

				viewervr.data.set_mesh_with_floor(V, F);
				viewervr.data.compute_normals();

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
                if(!added_stroke->has_points_on_mesh) {
                    hand_has_moved = false;
                    return true;
                }
                added_stroke->snap_to_vertices(vertex_boundary_markers);
                stroke_collection.push_back(*added_stroke);
                draw_all_strokes(viewervr);
            
		}
		else if (prev_tool_mode == REMOVE && stroke_was_removed) {

		}
		else if (prev_tool_mode == PULL && handleID != -1 && hand_has_moved) { //TODO: take care of hand_has_moved logic
		}

		prev_tool_mode = NONE;
		return true;
	}



	return true;
}

bool callback_mouse_down(ViewerVR& viewervr, int button, int modifier) {
	/*if(button == (int)ViewerVR::MouseButton::Right) {
		return false;
	}
	mouse_is_down = true;

	down_mouse_x = viewervr.current_mouse_x;
	down_mouse_y = viewervr.current_mouse_y;

	if(tool_mode == DRAW) { //Creating the first curve/mesh
		viewervr.data.clear();
		stroke_collection.clear();
		next_added_stroke_ID = 2;
		initial_stroke->strokeReset();
		initial_stroke->strokeAddSegment(down_mouse_x, down_mouse_y);
		skip_standardcallback = true;
	}*/ /*
	else if(tool_mode == ADD) { //Adding a new control curve onto an existing mesh
		added_stroke = new Stroke(V, F, viewer, next_added_stroke_ID);
		next_added_stroke_ID++;
		added_stroke->strokeAddSegmentAdd(down_mouse_x, down_mouse_y); //If the user starts outside of the mesh, consider the movement as navigation
		skip_standardcallback = true;
	} 
	else if(tool_mode == REMOVE) {
		double closest_dist = INFINITY;
		double current_closest = closest_dist;
		int tmp_handleID, closest_stroke_idx;
		handleID = -1;
		for(int i = 0; i < stroke_collection.size(); i++) {
			tmp_handleID = stroke_collection[i].selectClosestVertex(down_mouse_x, down_mouse_y, closest_dist);
			if((closest_dist < current_closest) && (tmp_handleID != -1)) {
				current_closest = closest_dist;
				handleID = tmp_handleID;
				closest_stroke_ID = stroke_collection[i].get_ID();
				closest_stroke_idx = i;
			}
		}

		if(handleID == -1) {//User clicked too far from any of the stroke vertices
			return false;
		}
		if(closest_stroke_ID == prev_closest_stroke_ID) {
			remove_stroke_clicked++;
		} else {
			remove_stroke_clicked = 1; //Start from 1
			prev_closest_stroke_ID = closest_stroke_ID;
		}

		//Redraw the original stroke and all added strokes, where the selected stroke is drawn in black.
		Eigen::MatrixXd init_points = initial_stroke->get3DPoints();
		viewer.data.set_points(init_points.topRows(init_points.rows() - 1), Eigen::RowVector3d(1, 0, 0));
		Eigen::MatrixXd added_points = stroke_collection[closest_stroke_idx].get3DPoints();
		viewer.data.add_points(added_points.topRows(added_points.rows() - 1), Eigen::RowVector3d(0, 0, 0));
		for(int i = 0; i < stroke_collection.size(); i++) {
			if(stroke_collection[i].get_ID() == closest_stroke_ID) {
				continue;
			}
			added_points = stroke_collection[i].get3DPoints();
			viewer.data.add_points(added_points.topRows(added_points.rows() - 1), stroke_collection[i].stroke_color);
		}


		if(remove_stroke_clicked == 2) { //Mechanism to force the user to click twice on the same stroke before removing it (safeguard)
			stroke_was_removed = true;
			stroke_collection[closest_stroke_idx].undo_stroke_add(vertex_boundary_markers); //Sets the vertex_boundary_markers for the vertices of this stroke to 0 again
			stroke_collection.erase(stroke_collection.begin() + closest_stroke_idx);
			remove_stroke_clicked = 0; //Reset
		}

		skip_standardcallback = true;
	} 
	else if(tool_mode == PULL) { //Dragging an existing curve
		select_dragging_handle(down_mouse_x, down_mouse_y);

		if(handleID == -1) {//User clicked too far from any of the stroke vertices
			return false;
		}
		if(closest_stroke_ID == -1) {
			CurveDeformation::startPullCurve(*initial_stroke, handleID, V.rows(), part_of_original_stroke);
		} else {
			CurveDeformation::startPullCurve(stroke_collection[closest_stroke_ID], handleID, V.rows(), part_of_original_stroke);
		}
		skip_standardcallback = true;
	} 
	else if(tool_mode == NAVIGATE) { //Navigate through the screen
		skip_standardcallback = false; //We do want to use the navigation functionality
	} 
	else if(tool_mode == CUT) {
		if(cut_stroke_already_drawn) { //clicked while cut stroke already drawn
			return true;
		}
		//clicked with no cut stroke drawn yet
		added_stroke = new Stroke(V, F, viewer, next_added_stroke_ID);
		next_added_stroke_ID++;
		added_stroke->strokeAddSegmentCut(down_mouse_x, down_mouse_y);
		skip_standardcallback = true;
	} 
	else if(tool_mode == EXTRUDE) {
		if(extrusion_base_already_drawn) { //clicked while the extrude base was already drawn
			added_stroke = new Stroke(V, F, viewer, next_added_stroke_ID);
			next_added_stroke_ID++;
			added_stroke->strokeAddSegmentExtrusionSilhouette(down_mouse_x, down_mouse_y);
		} else { //clicked with no extrude base yet
			extrusion_base = new Stroke(V, F, viewer, next_added_stroke_ID);
			next_added_stroke_ID++;
			extrusion_base->strokeAddSegmentExtrusionBase(down_mouse_x, down_mouse_y);
		}
		skip_standardcallback = true;
	}*/

	return skip_standardcallback; //Will make sure that we use standard navigation responses if we didn't do special actions and vice versa
}

bool callback_mouse_move(ViewerVR& viewervr, int mouse_x, int mouse_y) {
	if(!skip_standardcallback) {
		return false;
	}

	/*if(viewervr.down) { //Only consider it to be moving if the button was held down
		hand_has_moved = true;
	}


	if(tool_mode == DRAW && viewervr.down) { //If we're still holding the mouse down
		initial_stroke->strokeAddSegment(mouse_x, mouse_y);
		return true;
	}*/ /*
	else if(tool_mode == ADD && viewer.down) {
		last_add_on_mesh = added_stroke->strokeAddSegmentAdd(mouse_x, mouse_y);
		return true;
	} 
	else if(tool_mode == EXTRUDE && viewer.down) {
		if(extrusion_base_already_drawn) {
			added_stroke->strokeAddSegmentExtrusionSilhouette(mouse_x, mouse_y);
		} else {
			extrusion_base->strokeAddSegmentExtrusionBase(mouse_x, mouse_y);
		}
		return true;
	} 
	else if(tool_mode == PULL && viewer.down && handleID != -1) {
		double x = mouse_x;
		double y = viewer.core.viewport(3) - mouse_y;

		Eigen::Matrix4f modelview = viewer.core.view * viewer.core.model;
		int global_handleID;
		if(closest_stroke_ID == -1) {
			global_handleID = initial_stroke->get_vertex_idx_for_point(handleID);
		} else {
			global_handleID = stroke_collection[closest_stroke_ID].get_vertex_idx_for_point(handleID);
		}
		Eigen::RowVector3f pt1(viewer.data.V(global_handleID, 0), viewer.data.V(global_handleID, 1), viewer.data.V(global_handleID, 2));
		Eigen::RowVector3f pr;
		igl::project(pt1, modelview, viewer.core.proj, viewer.core.viewport, pr);
		Eigen::RowVector3d pt = igl::unproject(Eigen::Vector3f(x, y, pr[2]), modelview, viewer.core.proj, viewer.core.viewport).transpose().cast<double>();

		if(turnNr == 0) { //increase the number to smooth less often
			CurveDeformation::pullCurve(pt, V);
			if(dirty_boundary) { //Smooth an extra time if the boundary is dirty, because smoothing once with a dirty boundary results in a flat mesh
				for(int i = 0; i < 2; i++) {
						SurfaceSmoothing::smooth(V, F, vertex_boundary_markers, part_of_original_stroke, new_mapped_indices, sharp_edge, dirty_boundary);
				}
			}
				SurfaceSmoothing::smooth(V, F, vertex_boundary_markers, part_of_original_stroke, new_mapped_indices, sharp_edge, dirty_boundary);

			turnNr++;
		} else {
			turnNr++;
			if(turnNr == 4) {
				turnNr = 0;
			}
		}

		initial_stroke->update_Positions(V);
		for(int i = 0; i < stroke_collection.size(); i++) {
			stroke_collection[i].update_Positions(V);
		}

		viewer.data.set_mesh(V, F);
		viewer.data.compute_normals();
		draw_all_strokes(viewer);

		return true;
	} 
	else if(tool_mode == CUT && viewer.down) {
		if(!cut_stroke_already_drawn) {
			added_stroke->strokeAddSegmentCut(mouse_x, mouse_y);
		}
		return true;
	}*/

	return false;
}

bool callback_mouse_up(ViewerVR& viewervr, int button, int modifier) {
	if(!mouse_is_down) {
		return true;
	}
	mouse_is_down = false;

	/*if(tool_mode == DRAW) {
		if(initial_stroke->toLoop()) {//Returns false if the stroke only consists of 1 point (user just clicked)
			//Give some time to show the stroke
#ifdef _WIN32
			Sleep(200);
#else
			usleep(200000); 
#endif
			backside_vertex_map = initial_stroke->generate3DMeshFromStroke(vertex_boundary_markers, part_of_original_stroke);
			F = viewervr.data.F;
			V = viewervr.data.V;

			Eigen::MatrixXi EV, FE, EF;
			igl::edge_topology(V, F, EV, FE, EF);
			sharp_edge.resize(EV.rows());
			sharp_edge.setZero(); //Set all edges to smooth after initial draw

			dirty_boundary = true;

			for(int i = 0; i < initial_smooth_iter; i++) {
				SurfaceSmoothing::smooth(V, F, vertex_boundary_markers, part_of_original_stroke, new_mapped_indices, sharp_edge, dirty_boundary);
			}


			viewervr.data.set_mesh(V, F);
			viewervr.data.compute_normals();

			//Overlay the drawn stroke
			int strokeSize = (vertex_boundary_markers.array() > 0).count();
			Eigen::MatrixXd strokePoints = V.block(0, 0, strokeSize, 3);
			viewervr.data.set_points(strokePoints, Eigen::RowVector3d(1, 0, 0)); //Displays dots
			viewervr.data.set_stroke_points(igl::cat(1, strokePoints, (Eigen::MatrixXd) V.row(0)));

		}
		skip_standardcallback = false;
	}*//*
	else if(tool_mode == ADD) {
		dirty_boundary = true;
		if(!added_stroke->has_points_on_mesh) {
			hand_has_moved = false;
			return true;
		}
		added_stroke->snap_to_vertices(vertex_boundary_markers);
        stroke_collection.push_back(*added_stroke);
		draw_all_strokes(viewer);
	} 
	else if(tool_mode == REMOVE && stroke_was_removed) { //Only redraw if we actually removed a stroke (otherwise we draw unnecessary)
		stroke_was_removed = false; //Reset
		dirty_boundary = true;

		draw_all_strokes(viewer);
	} 
	else if(tool_mode == PULL && handleID != -1 && hand_has_moved) {
		for(int i = 0; i < 2; i++) {
			     SurfaceSmoothing::smooth(V, F, vertex_boundary_markers, part_of_original_stroke, new_mapped_indices, sharp_edge, dirty_boundary);
		}

		for(int i = 0; i < stroke_collection.size(); i++) {
			stroke_collection[i].update_Positions(V);
		}

		viewer.data.set_mesh(V, F);
		viewer.data.compute_normals(); //TODO: might need to use igl::per_face_normals(V, F, N_Faces); viewer.data.set_normals(N_Faces);

		draw_all_strokes(viewer);

	} 
	else if(tool_mode == CUT) {
		if(!added_stroke->has_points_on_mesh) {
			hand_has_moved = false;
			return true;
		}
		if(cut_stroke_already_drawn) { //User had already drawn the cut stroke and has now drawn the final stroke for removing the part
			dirty_boundary = true;
			added_stroke->append_final_point();
			added_stroke->toLoop();
			MeshCut::cut(V, F, vertex_boundary_markers, part_of_original_stroke, new_mapped_indices, sharp_edge, *added_stroke);
			stroke_collection.push_back(*added_stroke);

			initial_stroke->update_vert_bindings(new_mapped_indices, vertex_boundary_markers);//Don't test if the initial one dies, cause then we have mayhem anyway? TODO
			
			int nr_removed = 0, original_collection_size = stroke_collection.size();
			for(int i = 0; i < original_collection_size; i++) {
				if(!stroke_collection[i - nr_removed].update_vert_bindings(new_mapped_indices, vertex_boundary_markers)) {
					//Stroke dies, don't need to do stroke.undo_stroke_add, cause all its vertices also cease to exist
					stroke_collection.erase(stroke_collection.begin() + i - nr_removed);
					nr_removed++;
					dirty_boundary = true;
					continue; //Go to the next stroke, don't update this ones' positions
				}
			}

			for(int i = 0; i < 2; i++) {
				SurfaceSmoothing::smooth(V, F, vertex_boundary_markers, part_of_original_stroke, new_mapped_indices, sharp_edge, dirty_boundary);
			}

			//Update the stroke positions after smoothing, in case their positions have changed (although they really shouldn't)
			initial_stroke->update_Positions(V);
			for(int i = 0; i < stroke_collection.size(); i++) {
				stroke_collection[i].update_Positions(V);
			}

			viewer.data.clear();
			viewer.data.set_mesh(V, F);
			igl::per_face_normals(V, F, N_Faces);
			viewer.data.set_normals(N_Faces);
			viewer.core.align_camera_center(V, F);

			cut_stroke_already_drawn = false; //Reset
			draw_all_strokes(viewer);
		} else { //We're finished drawing the cut stroke, prepare for when user draws the final stroke to remove the part
			cut_stroke_already_drawn = true;
		}
	} 
	else if(tool_mode == EXTRUDE) {
		if(extrusion_base_already_drawn) { //User has drawn the silhouette stroke for extrusion
			dirty_boundary = true;
			cout << "mouse released after extrusion silhouette drawn" << endl;
			added_stroke->toLoop();
			MeshExtrusion::extrude_main(V, F, vertex_boundary_markers, part_of_original_stroke, new_mapped_indices, sharp_edge, base_surface_path, *added_stroke, *extrusion_base, base_model, base_view, base_proj, base_viewport);
			stroke_collection.push_back(*extrusion_base);
			stroke_collection.push_back(*added_stroke);

			initial_stroke->update_vert_bindings(new_mapped_indices, vertex_boundary_markers);

			int nr_removed = 0, original_collection_size = stroke_collection.size();
			for(int i = 0; i < original_collection_size - 1; i++) { //Skip the newly added stroke since it is already updated inside extrude_main()
				if(!stroke_collection[i - nr_removed].update_vert_bindings(new_mapped_indices, vertex_boundary_markers)) {
					//Stroke dies, don't need to do stroke.undo_stroke_add, cause all its vertices also cease to exist (or in case of non-loop strokes that have a middle portion removed, the undo_stroke_add is done inside of the update_vert_bindings
					stroke_collection.erase(stroke_collection.begin() + i - nr_removed);
					nr_removed++;
					dirty_boundary = true;
					continue; //Go to the next stroke, don't update this ones' positions
				}
			}

			for(int i = 0; i < 3; i++) {
				SurfaceSmoothing::smooth(V, F, vertex_boundary_markers, part_of_original_stroke, new_mapped_indices, sharp_edge, dirty_boundary);
			}

			//Update the stroke positions after smoothing, in case their positions have changed (although they really shouldn't)
			initial_stroke->update_Positions(V);
			for(int i = 0; i < stroke_collection.size()-1; i++) {
				stroke_collection[i].update_Positions(V);
			}

			viewer.data.clear();
			viewer.data.set_mesh(V, F);

			igl::per_face_normals(V, F, N_Faces);
			viewer.data.set_normals(N_Faces);
			viewer.core.align_camera_center(V, F);

			extrusion_base_already_drawn = false; //Reset
			draw_all_strokes(viewer);
		} else { //mouse released after extrusion base drawn
			if(!extrusion_base->has_points_on_mesh) {
				hand_has_moved = false;
				return true;
			}
			//extrusion_base->resample_all(); //This will shrink the drawn stroke. Might result in no face being contained inside the stroke
			dirty_boundary = true;
			extrusion_base->toLoop();
			extrusion_base_already_drawn = true;
			MeshExtrusion::extrude_prepare(*extrusion_base, base_surface_path); //Don't need to update all strokes here, since it didn't remove any vertices

			base_model = extrusion_base->viewer.core.model;
			base_view = extrusion_base->viewer.core.view;
			base_proj = extrusion_base->viewer.core.proj;
			base_viewport = extrusion_base->viewer.core.viewport;

			draw_all_strokes(viewer);
			draw_extrusion_base(viewer); //Need to draw the extrusion base separately, since it isn't added to the stroke_collection yet.
		}
	}*/

	hand_has_moved = false;
	return skip_standardcallback;
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
	ViewerVR viewervr;
	igl::read_triangle_mesh("../data/cube.off", V, F);
	//cout << V << endl << endl << F << endl;
	viewervr.data.set_mesh_with_floor(V, F);
/*	viewervr.callback_key_down = callback_key_down;
	viewervr.callback_mouse_down = callback_mouse_down;
	viewervr.callback_mouse_move = callback_mouse_move;
	viewervr.callback_mouse_up = callback_mouse_up;
	viewervr.core.point_size = 15;*/
	//viewer.callback_load_mesh = callback_load_mesh;

	/*viewer.callback_init = [&](igl::viewer::Viewer& viewer) {
		// Add new group
		viewer.ngui->addGroup("Inflation");

		// Expose a variable directly ...
		viewer.ngui->addVariable("Vertex Weights", SurfaceSmoothing::vertex_weight);
		viewer.ngui->addVariable("Edge Weights", SurfaceSmoothing::edge_weight);


		// Expose a variable directly ...
		viewer.ngui->addVariable("Initial smoothing iterations", initial_smooth_iter);


		// Add a button
		viewer.ngui->addButton("Perform 1 smoothing iteration", [&viewer]() {
			SurfaceSmoothing::smooth(V, F, vertex_boundary_markers, part_of_original_stroke, new_mapped_indices, sharp_edge, dirty_boundary);
			viewer.data.set_mesh(V, F);
			viewer.data.compute_normals();


			for(int i = 0; i < stroke_collection.size(); i++) {
				stroke_collection[i].update_Positions(V);
			}

			//Overlay the updated stroke
			initial_stroke->update_Positions(V);
			Eigen::MatrixXd init_points = initial_stroke->get3DPoints();
			viewer.data.set_points(init_points.topRows(init_points.rows() - 1), Eigen::RowVector3d(1, 0, 0));
			viewer.data.set_stroke_points(init_points);

			viewer.data.set_edges(Eigen::MatrixXd(), Eigen::MatrixXi(), Eigen::RowVector3d(0, 0, 1));
			Eigen::MatrixXd added_points;
			int points_to_hold_back;
			for(int i = 0; i < stroke_collection.size(); i++) {
				added_points = stroke_collection[i].get3DPoints();
				points_to_hold_back = 1 + !stroke_collection[i].is_loop;
				viewer.data.add_points(added_points.topRows(added_points.rows() - 1), stroke_collection[i].stroke_color);
				viewer.data.add_edges(added_points.block(0, 0, added_points.rows() - points_to_hold_back, 3), added_points.block(1, 0, added_points.rows() - points_to_hold_back, 3), stroke_collection[i].stroke_color);// Eigen::RowVector3d(0, 0, 1));
			}

		});

		viewer.ngui->addGroup("Curve deformation");
		viewer.ngui->addVariable<bool>("Smooth deformation", CurveDeformation::smooth_deform_mode);

		// call to generate menu
		viewer.screen->performLayout();
		return false;
	};*/

	//Init stroke selector
	initial_stroke = new Stroke(V, F, viewervr, 0);
	if(argc == 2) {
		// Read mesh
	//	igl::readOFF(argv[1], V, F);
		//	callback_load_mesh(viewer, argv[1]);
	} else {
		// Read mesh
		//callback_load_mesh(viewer, "../data/cube.off");
	}

//	callback_key_down(viewervr, '1', 0);

	//viewer.core.align_camera_center(V);
    viewervr.init();
	viewervr.callback_button_down = button_down;

	viewervr.launch();
}

