#include "Stroke.h"
#include <igl/unproject_onto_mesh.h>
#include <igl/unproject.h>
#include <igl/unproject_ray.h>
#include <igl/triangle/triangulate.h>
#include <algorithm> 
#include <igl/per_face_normals.h>
#include <igl/per_vertex_normals.h>
#include <igl/dijkstra.h>
#include <igl/edge_topology.h>
#include "Plane.h"
using namespace igl;
using namespace std;

std::chrono::steady_clock::time_point _time2, _time1;
//Eigen::MatrixXi EV, FE, EF;

Stroke::Stroke(const Eigen::MatrixXd &V_, const Eigen::MatrixXi &F_, igl::viewer::Viewer &v, int stroke_ID_) :
	V(V_),
	F(F_),
	viewer(v), 
	stroke_ID(stroke_ID_) {
	stroke2DPoints = Eigen::MatrixX2d::Zero(1, 2);
	stroke3DPoints = Eigen::MatrixX3d::Zero(1, 3);
	stroke_edges = Eigen::MatrixXi::Zero(0, 2);
	_time1 = std::chrono::high_resolution_clock::now();
	closest_vert_bindings.clear();
	has_points_on_mesh = false;
	stroke_color = Eigen::RowVector3d(0.8*(rand() / (double)RAND_MAX), 0.8*(rand() / (double)RAND_MAX), 0.8*(rand() / (double)RAND_MAX));
}

Stroke::Stroke(const Stroke& origin) :
	V(origin.V),
	F(origin.F),
	viewer(origin.viewer),
	stroke_ID(origin.stroke_ID),
	stroke2DPoints(origin.stroke2DPoints),
	stroke3DPoints(origin.stroke3DPoints),
	stroke_edges(origin.stroke_edges),
	closest_vert_bindings(origin.closest_vert_bindings),
	has_points_on_mesh(origin.has_points_on_mesh),
	stroke_color(origin.stroke_color),
	dep(origin.dep),
	added_stroke_final_vertices(origin.added_stroke_final_vertices),
	is_loop(origin.is_loop)
	{
	_time1 = std::chrono::high_resolution_clock::now();
}

Stroke& Stroke::operator=(Stroke other) {
	swap(other);
	return *this;
}

void Stroke::swap(Stroke & tmp) {//The pointers to V and F will always be the same for all stroke instances, so no need to copy
	std::swap(this->viewer, tmp.viewer);
	std::swap(this->stroke_ID, tmp.stroke_ID);
	std::swap(this->stroke2DPoints, tmp.stroke2DPoints);
	std::swap(this->stroke3DPoints, tmp.stroke3DPoints);
	std::swap(this->stroke_edges, tmp.stroke_edges);
	std::swap(this->closest_vert_bindings, tmp.closest_vert_bindings);
	std::swap(this->has_points_on_mesh, tmp.has_points_on_mesh);
	std::swap(this->stroke_color, tmp.stroke_color);
	std::swap(this->dep, tmp.dep);
	std::swap(this->added_stroke_final_vertices, tmp.added_stroke_final_vertices);
	std::swap(this->is_loop, tmp.is_loop);
}

Stroke::~Stroke() {}

//To be used exclusively for drawing strokes for INITIAL MESHES (no extrusion)
void Stroke::strokeAddSegment(int mouse_x, int mouse_y) {
	//OpenGL has origin at left bottom, window(s) has origin at left top
	double x = mouse_x;
	double y = viewer.core.viewport(3) - mouse_y;
	if (!empty2D() && x == stroke2DPoints(stroke2DPoints.rows() - 1, 0) && y == stroke2DPoints(stroke2DPoints.rows() - 1, 1)) { //Check that the point is new compared to last time
		return;
	}

	if (!empty2D()) {
		_time2 = std::chrono::high_resolution_clock::now();
		auto timePast = std::chrono::duration_cast<std::chrono::nanoseconds>(_time2 - _time1).count();
        if (timePast < 10000000) { //Don't add another segment before x nanoseconds
			return;
		}
	}

	Eigen::Matrix4f modelview = viewer.core.view * viewer.core.model;
	Eigen::RowVector3d pt(0, 0, 0);

	pt = igl::unproject(Eigen::Vector3f(x, y, 0.0f), modelview, viewer.core.proj, viewer.core.viewport).transpose().cast<double>();
	if (stroke2DPoints.rows() == 1 && empty2D()) { //Add first point
		stroke2DPoints.row(0) << x, y;
		stroke3DPoints.row(0) << pt[0], pt[1], pt[2];
	} else {
		stroke2DPoints.conservativeResize(stroke2DPoints.rows() + 1, stroke2DPoints.cols());
		stroke2DPoints.row(stroke2DPoints.rows() - 1) << x, y;

		stroke3DPoints.conservativeResize(stroke3DPoints.rows() + 1, stroke3DPoints.cols());
		stroke3DPoints.row(stroke3DPoints.rows() - 1) << pt[0], pt[1], pt[2];

		stroke_edges.conservativeResize(stroke_edges.rows() + 1, stroke_edges.cols());
		stroke_edges.row(stroke_edges.rows() - 1) << stroke2DPoints.rows() - 2, stroke2DPoints.rows() - 1;
	}
	closest_vert_bindings.push_back(stroke3DPoints.rows() - 1);
	//using set_stroke_points will remove all previous strokes, using add_stroke_points might create duplicates
	viewer.data.set_stroke_points(stroke3DPoints);

	_time1 = std::chrono::high_resolution_clock::now(); //restart the "start" timer
}

bool Stroke::strokeAddSegmentAdd(int mouse_x, int mouse_y) {
	bool result = false;
	//OpenGL has origin at left bottom, window(s) has origin at left top
	double x = mouse_x;
	double y = viewer.core.viewport(3) - mouse_y;
	if(!empty2D() && x == stroke2DPoints(stroke2DPoints.rows() - 1, 0) && y == stroke2DPoints(stroke2DPoints.rows() - 1, 1)) { //Check that the point is new compared to last time
		return result;
	}

	if(!empty2D()) {
		_time2 = std::chrono::high_resolution_clock::now();
		auto timePast = std::chrono::duration_cast<std::chrono::nanoseconds>(_time2 - _time1).count();
		if(timePast < 10000000) { //Don't add another segment before x nanoseconds
			return result;
		}
	}

	Eigen::Matrix4f modelview = viewer.core.view * viewer.core.model;
	Eigen::RowVector3d pt(0, 0, 0);
	int faceID = -1;

	Eigen::Vector3f bc;
	if(igl::unproject_onto_mesh(Eigen::Vector2f(x, y), modelview, viewer.core.proj, viewer.core.viewport, V, F, faceID, bc)) {
		pt = V.row(F(faceID, 0))*bc(0) + V.row(F(faceID, 1))*bc(1) + V.row(F(faceID, 2))*bc(2);
		double cur_dist, min_dist = INFINITY;
		int closest_vert_idx;
		for(int i = 0; i < 3; i++) {
			cur_dist = (pt - V.row(F(faceID, i))).norm();
			if(cur_dist < min_dist) {
				min_dist = cur_dist;
				closest_vert_idx = F(faceID, i);
			}
		}

		if(closest_vert_bindings.size()==0 || closest_vert_idx != closest_vert_bindings.back()) {
			closest_vert_bindings.push_back(closest_vert_idx);

		}

		has_points_on_mesh = true;

		//Only add point when it's on the mesh
		if(stroke2DPoints.rows() == 1 && empty2D()) { //Add first point
			stroke2DPoints.row(0) << x, y;
			stroke3DPoints.row(0) << pt[0], pt[1], pt[2];
		} else {
			stroke2DPoints.conservativeResize(stroke2DPoints.rows() + 1, stroke2DPoints.cols());
			stroke2DPoints.row(stroke2DPoints.rows() - 1) << x, y;

			stroke3DPoints.conservativeResize(stroke3DPoints.rows() + 1, stroke3DPoints.cols());
			stroke3DPoints.row(stroke3DPoints.rows() - 1) << pt[0], pt[1], pt[2];

			stroke_edges.conservativeResize(stroke_edges.rows() + 1, stroke_edges.cols());
			stroke_edges.row(stroke_edges.rows() - 1) << stroke2DPoints.rows() - 2, stroke2DPoints.rows() - 1;
		}
		//using set_stroke_points will remove all previous strokes, using add_stroke_points might create duplicates
		viewer.data.add_edges(stroke3DPoints.block(0, 0, stroke3DPoints.rows() - 1, 3), stroke3DPoints.block(1, 0, stroke3DPoints.rows() - 1, 3), Eigen::RowVector3d(0, 1, 0));
		result = true;
	}

	_time1 = std::chrono::high_resolution_clock::now(); //restart the "start" timer
	return result;
}

bool Stroke::strokeAddSegmentCut(int mouse_x, int mouse_y) {
	bool result = false;
	//OpenGL has origin at left bottom, window(s) has origin at left top
	double x = mouse_x;
	double y = viewer.core.viewport(3) - mouse_y;
	if(!empty2D() && x == stroke2DPoints(stroke2DPoints.rows() - 1, 0) && y == stroke2DPoints(stroke2DPoints.rows() - 1, 1)) { //Check that the point is new compared to last time
		return result;
	}

	if(!empty2D()) {
		_time2 = std::chrono::high_resolution_clock::now();
		auto timePast = std::chrono::duration_cast<std::chrono::nanoseconds>(_time2 - _time1).count();
		if(timePast < 10000000) { //Don't add another segment before x nanoseconds
			return result;
		}
	}

	Eigen::Matrix4f modelview = viewer.core.view * viewer.core.model;
	Eigen::RowVector3d pt(0, 0, 0);
	int faceID = -1;

	Eigen::Vector3f bc;
	if(igl::unproject_onto_mesh(Eigen::Vector2f(x, y), modelview, viewer.core.proj, viewer.core.viewport, V, F, faceID, bc)) {
		pt = V.row(F(faceID, 0))*bc(0) + V.row(F(faceID, 1))*bc(1) + V.row(F(faceID, 2))*bc(2);
		has_points_on_mesh = true;

		//Only add point when it's on the mesh
		if(stroke2DPoints.rows() == 1 && empty2D()) { //Add first point
			stroke2DPoints.row(0) << x, y;
			stroke3DPoints.row(0) << pt[0], pt[1], pt[2];
		} else {
			stroke2DPoints.conservativeResize(stroke2DPoints.rows() + 1, stroke2DPoints.cols());
			stroke2DPoints.row(stroke2DPoints.rows() - 1) << x, y;

			stroke3DPoints.conservativeResize(stroke3DPoints.rows() + 1, stroke3DPoints.cols());
			stroke3DPoints.row(stroke3DPoints.rows() - 1) << pt[0], pt[1], pt[2];

			stroke_edges.conservativeResize(stroke_edges.rows() + 1, stroke_edges.cols());
			stroke_edges.row(stroke_edges.rows() - 1) << stroke2DPoints.rows() - 2, stroke2DPoints.rows() - 1;
		}
		//using set_stroke_points will remove all previous strokes, using add_stroke_points might create duplicates
		viewer.data.add_points(stroke3DPoints.row(stroke3DPoints.rows() - 1), Eigen::RowVector3d(1, 0, 1));
	//	viewer.data.add_edges(stroke3DPoints.block(0, 0, stroke3DPoints.rows() - 1, 3), stroke3DPoints.block(1, 0, stroke3DPoints.rows() - 1, 3), Eigen::RowVector3d(0, 1, 0));
		result = true;
		just_came_from_mesh = true;
	} else if(stroke2DPoints.rows()>1 && just_came_from_mesh) { //We need to add the final point of the stroke, even though it is off the mesh (needed in order to wrap around to the backside)
		cut_stroke_final_point = igl::unproject(Eigen::Vector3f(x, y, 0.95*dep), modelview, viewer.core.proj, viewer.core.viewport).transpose().cast<double>();
		cut_stroke_final_point_2D = Eigen::RowVector2d(x, y);
		just_came_from_mesh = false;
	}

	_time1 = std::chrono::high_resolution_clock::now(); //restart the "start" timer
	return result;
}

void Stroke::append_final_point() {
	stroke2DPoints.conservativeResize(stroke2DPoints.rows() + 1, stroke2DPoints.cols());
	stroke2DPoints.row(stroke2DPoints.rows() - 1) << cut_stroke_final_point_2D[0], cut_stroke_final_point_2D[1];

	stroke3DPoints.conservativeResize(stroke3DPoints.rows() + 1, stroke3DPoints.cols());
	stroke3DPoints.row(stroke3DPoints.rows() - 1) << cut_stroke_final_point[0], cut_stroke_final_point[1], cut_stroke_final_point[2];

	stroke_edges.conservativeResize(stroke_edges.rows() + 1, stroke_edges.cols());
	stroke_edges.row(stroke_edges.rows() - 1) << stroke2DPoints.rows() - 2, stroke2DPoints.rows() - 1;

	viewer.data.add_edges(stroke3DPoints.block(0, 0, stroke3DPoints.rows() - 1, 3), stroke3DPoints.block(1, 0, stroke3DPoints.rows() - 1, 3), Eigen::RowVector3d(0, 1, 0));
}


//For drawing extrusion strokes. Need to start on the existing mesh
void Stroke::strokeAddSegmentExtrusion(int mouse_x, int mouse_y) {
	//OpenGL has origin at left bottom, window(s) has origin at left top
	double x = mouse_x;
	double y = viewer.core.viewport(3) - mouse_y;
	if (!empty2D() && x == stroke2DPoints(stroke2DPoints.rows() - 1, 0) && y == stroke2DPoints(stroke2DPoints.rows() - 1, 1)) { //Check that the point is new compared to last time
		return;
	}

	if (!empty2D()) {
		_time2 = std::chrono::high_resolution_clock::now();
		auto timePast = std::chrono::duration_cast<std::chrono::nanoseconds>(_time2 - _time1).count();
		if (timePast < 10000000) { //Don't add another segment before x nanoseconds
			return;
		}
	}

	Eigen::Matrix4f modelview = viewer.core.view * viewer.core.model;
	Eigen::RowVector3d pt(0, 0, 0);
	int faceID = -1;

	if (dep < 0) { //first time
		Eigen::Vector3f bc;
		if (igl::unproject_onto_mesh(Eigen::Vector2f(x, y), modelview, viewer.core.proj, viewer.core.viewport, V, F, faceID, bc)) {
			pt = V.row(F(faceID, 0))*bc(0) + V.row(F(faceID, 1))*bc(1) + V.row(F(faceID, 2))*bc(2);
			Eigen::Vector3f proj = igl::project(pt.transpose().cast<float>().eval(), modelview, viewer.core.proj, viewer.core.viewport);
			dep = proj[2]; //0.85 is at center of meshes?
		}
	}
	pt = igl::unproject(Eigen::Vector3f(x, y, 0.95*dep), modelview, viewer.core.proj, viewer.core.viewport).transpose().cast<double>();

	if (stroke2DPoints.rows() == 1 && empty2D()) { //Add first point
		stroke2DPoints.row(0) << x, y;
		stroke3DPoints.row(0) << pt[0], pt[1], pt[2];
	} else {
		stroke2DPoints.conservativeResize(stroke2DPoints.rows() + 1, stroke2DPoints.cols());
		stroke2DPoints.row(stroke2DPoints.rows() - 1) << x, y;

		stroke3DPoints.conservativeResize(stroke3DPoints.rows() + 1, stroke3DPoints.cols());
		stroke3DPoints.row(stroke3DPoints.rows() - 1) << pt[0], pt[1], pt[2];

		stroke_edges.conservativeResize(stroke_edges.rows() + 1, stroke_edges.cols());
		stroke_edges.row(stroke_edges.rows() - 1) << stroke2DPoints.rows() - 2, stroke2DPoints.rows() - 1;
	}
	//using set_stroke_points will remove all previous strokes, using add_stroke_points might create duplicates
	viewer.data.set_stroke_points(stroke3DPoints);

	_time1 = std::chrono::high_resolution_clock::now(); //restart the "start" timer
}

void Stroke::strokeReset() {
	stroke2DPoints.resize(1, 2);
	stroke2DPoints.setZero();
	stroke3DPoints.resize(1, 3);
	stroke3DPoints.setZero();
	stroke_edges.resize(0, 2);
	stroke_edges.setZero();
	dep = -1;
	_time1 = std::chrono::high_resolution_clock::now();
	closest_vert_bindings.clear();
}

bool Stroke::toLoop() {
	if(stroke2DPoints.rows() > 2) { //Don't do anything if we have only 1 line segment
		stroke_edges.conservativeResize(stroke_edges.rows() + 1, stroke_edges.cols());
		stroke_edges.row(stroke_edges.rows() - 1) << stroke_edges.rows() - 1, 0; //Add an edge from the last vertex to the first
		//using set_stroke_points will remove all previous strokes, using add_stroke_points might create duplicates
        
        stroke3DPoints.conservativeResize(stroke3DPoints.rows()+1, stroke3DPoints.cols());
        stroke3DPoints.row(stroke3DPoints.rows()-1) << stroke3DPoints.row(0);
		closest_vert_bindings.push_back(0);
		viewer.data.set_stroke_points(stroke3DPoints);
		is_loop = true;
		return true;
	}
	is_loop = false;
	return false;
}

unordered_map<int, int> Stroke::generate3DMeshFromStroke(Eigen::VectorXi &vertex_boundary_markers, Eigen::VectorXi &part_of_original_stroke) {
	counter_clockwise(); //Ensure the stroke is counter-clockwise, handy later

	Eigen::MatrixX2d original_stroke2DPoints = stroke2DPoints;
	stroke2DPoints = resample_stroke(original_stroke2DPoints);

	Eigen::MatrixXd V2_tmp, V2;
	Eigen::MatrixXi F2, F2_back, vertex_markers, edge_markers;
	igl::triangle::triangulate((Eigen::MatrixXd) stroke2DPoints, stroke_edges, Eigen::MatrixXd(0,0), Eigen::MatrixXi::Constant(stroke2DPoints.rows(),1,1), Eigen::MatrixXi::Constant(stroke_edges.rows(),1,1), "Qq25", V2_tmp, F2, vertex_markers, edge_markers); //Capital Q silences triangle's output in cmd line. Also retrieves markers to indicate whether or not an edge/vertex is on the mesh boundary
	V2 = Eigen::MatrixXd::Zero(V2_tmp.rows(), V2_tmp.cols() + 1);
  
    //zero mean in x and y
    V2_tmp.col(0) = V2_tmp.col(0).array() - V2_tmp.col(0).mean();
    V2_tmp.col(1) = V2_tmp.col(1).array() - V2_tmp.col(1).mean();
    
	V2.block(0,0, V2_tmp.rows(), 2) = V2_tmp;


	generate_backfaces(F2, F2_back);
	int nr_boundary_vertices = (vertex_markers.array() == 1).count(); //check how many vertices were marked as 1 (boundary) TODO: CHECK THAT 1 IS INDEED BOUNDARY
	int nr_boundary_edges = (edge_markers.array() == 1).count(); //check how many edges were marked as 1 (boundary)
	
	int original_size = V2.rows();
	V2.conservativeResize(2*original_size - nr_boundary_vertices, V2.cols()); //Increase size, such that boundary vertices only get included one
	unordered_map<int, int> backside_vertex_map;
	int count = 0;
	for(int i = 0; i < original_size; i++) {
		if(vertex_markers(i) != 1) { //If the vertex is NOT on the boundary, duplicate it
			V2.row(original_size + count) << V2.row(i);
			backside_vertex_map.insert({i, original_size + count});
			backside_vertex_map.insert({original_size + count, i});
			count++;
		} else {
			backside_vertex_map.insert({i, i});
		}
	}
    
    V2.col(0) = V2.col(0).array() - V2.col(0).mean();
    V2.col(1) = V2.col(1).array() - V2.col(1).mean();

	//Add backside faces, using original vertices for boundary and copied vertices for inside
	//Duplicate all faces
	original_size = F2_back.rows();
	F2.conservativeResize(original_size*2, F2.cols());
	for(int i = 0; i < original_size; i++) {
		for(int j = 0; j < 3; j++) {
			auto got = backside_vertex_map.find(F2_back(i, j));
			F2(original_size + i, j) = got->second;
		}
	}

	Eigen::MatrixXd N_Faces, N_Vertices;
	igl::per_face_normals(V2, F2, N_Faces);
	igl::per_vertex_normals(V2, F2, PER_VERTEX_NORMALS_WEIGHTING_TYPE_UNIFORM, N_Vertices);
	
	vertex_boundary_markers.resize(V2.rows());
    part_of_original_stroke.resize(V2.rows());
	for(int i = 0; i < V2.rows(); i++) {
		if(i >= vertex_markers.rows()) { //vertex can't be boundary (it's on backside)
            V2.row(i) = V2.row(i) + 0.1*N_Vertices.row(i);
			vertex_boundary_markers[i] = 0;
            part_of_original_stroke[i] = 0;
		} else {
			if(vertex_markers(i) == 1) { //Don't change boundary vertices
				vertex_boundary_markers[i] = 1;
                part_of_original_stroke[i] = 1;
				stroke3DPoints.row(i) = V2.row(i);
				continue;
			}
            V2.row(i) = V2.row(i) + 0.1*N_Vertices.row(i);
			vertex_boundary_markers[i] = 0;
            part_of_original_stroke[i] = 0;
		}
	}

	viewer.data.clear();
	viewer.data.set_mesh(V2, F2);
	igl::per_face_normals(V2, F2, N_Faces);

	viewer.data.set_normals(N_Faces);
	viewer.core.align_camera_center(viewer.data.V);
	return backside_vertex_map;
}

void Stroke::counter_clockwise() {
	double total_area = 0;
	Eigen::Vector2d prev, next;
	prev = (Eigen::Vector2d) stroke2DPoints.row(stroke2DPoints.rows() - 1);
	for(int i = 0; i < stroke2DPoints.rows(); i++) {
		next = (Eigen::Vector2d) stroke2DPoints.row(i);
		total_area += (prev[1] + next[1]) * (next[0] - prev[0]);
		prev = next;
	}

	if(total_area > 0) { //reverse the vector
		stroke2DPoints = stroke2DPoints.colwise().reverse().eval();
	}

}

Eigen::MatrixX2d Stroke::resample_stroke(Eigen::MatrixX2d & original_stroke2DPoints) {
	Eigen::MatrixX2d new_stroke2DPoints = Eigen::MatrixX2d::Zero(original_stroke2DPoints.rows(), 2);
	int nr_iterations = max(5.0, original_stroke2DPoints.rows() / 4.0);
	for(int i = 0; i < nr_iterations; i++) {
		move_to_middle(original_stroke2DPoints, new_stroke2DPoints);
		move_to_middle(new_stroke2DPoints, original_stroke2DPoints); //TODO: is this to save us a copy? basically performing an extra move_to_middle step
	}
	return new_stroke2DPoints;
}

void Stroke::move_to_middle(Eigen::MatrixX2d &positions, Eigen::MatrixX2d &new_positions) {
	int n = positions.rows();
	//Do seperately for i=0, because modulo gives -1
	Eigen::Vector2d prev = positions.row(n - 1);
	Eigen::Vector2d cur = positions.row(0);
	Eigen::Vector2d next = positions.row(1);
	new_positions(0, 0) = (cur[0] * 2 + prev[0] + next[0]) / 4;
	new_positions(0, 1) = (cur[1] * 2 + prev[1] + next[1]) / 4;

	for(int i = 1; i < n; i++) {
		prev = positions.row((i - 1) % n);
		cur = positions.row(i%n);
		next = positions.row((i + 1) % n);

		new_positions(i, 0) = (cur[0] * 2 + prev[0] + next[0]) / 4;
		new_positions(i, 1) = (cur[1] * 2 + prev[1] + next[1]) / 4;
	}
}

/*double Stroke::total_stroke_length() {
	double total_length = 0;
	for(int i = 0; i < stroke2DPoints.rows(); i++) {
		total_length += (stroke2DPoints.row(i) - stroke2DPoints.row((i + 1) % stroke2DPoints.rows())).norm();
	}
	return total_length;
}*/

void Stroke::generate_backfaces(Eigen::MatrixXi &faces, Eigen::MatrixXi &back_faces) {
	back_faces = faces.rowwise().reverse().eval();
}

Eigen::MatrixX3d Stroke::get3DPoints() {
	return stroke3DPoints;
}

int Stroke::get_vertex_idx_for_point(int pt_idx) {
	return closest_vert_bindings[pt_idx];
}

vector<int> Stroke::get_closest_vert_bindings() {
	return closest_vert_bindings;
}

int Stroke::selectClosestVertex(int mouse_x, int mouse_y, double& closest_distance) {
	double x = mouse_x;
	double y = viewer.core.viewport(3) - mouse_y;

	Eigen::Matrix4f modelview = viewer.core.view * viewer.core.model;
	Eigen::RowVector3d pt, vert, clicked_pt(x,y,0);
	double closest_dist = INFINITY, dist;
	int closest_ID;
	for(int i = 0; i < stroke3DPoints.rows(); i++) { //ignore the last stroke points because it's a duplicate of the first
		vert = stroke3DPoints.row(i).transpose();
		igl::project(vert, modelview, viewer.core.proj, viewer.core.viewport, pt); //project the boundary vertex and store in pt
		dist = (pt - clicked_pt).squaredNorm();
		if(dist < closest_dist) {
			closest_dist = dist;
			closest_ID = i;
		}
	}

	closest_dist = sqrt(closest_dist);
	double stroke_diag = compute_stroke_diag();
	closest_distance = closest_dist;
	if(closest_dist > 0.08*stroke_diag) { //Don't do pulling when the user clicks too far away from any curve (treat as navigation movement)
		return -1;
	}
	return closest_ID;
}

double Stroke::compute_stroke_diag() {
	Eigen::Vector3d maxBB = stroke3DPoints.colwise().maxCoeff();
	Eigen::Vector3d minBB = stroke3DPoints.colwise().minCoeff();
	return (maxBB - minBB).norm();
}

void Stroke::update_Positions(Eigen::MatrixXd V) {
	for(int i = 0; i < stroke3DPoints.rows()-1; i++) {
		stroke3DPoints.row(i) = V.row(closest_vert_bindings[i]);
	}
	stroke3DPoints.row(stroke3DPoints.rows() - 1) = stroke3DPoints.row(0); //The last vertex in the stroke is the first point again (and not the first inside mesh vertex)
}

void Stroke::snap_to_vertices(Eigen::VectorXi &vertex_boundary_markers) {
	Eigen::VectorXd min_dist;
	Eigen::VectorXi previous;
	vector<int> result_path;
	int prev = -1;
	vector<vector<int>> adj_list;
	adjacency_list(F, adj_list);

	//Determine whether the stroke is a loop or not
	if((stroke3DPoints.row(0) - stroke3DPoints.row(stroke3DPoints.rows() - 1)).norm() < compute_stroke_diag() / 5.0) {
		is_loop = true;
	} else {
		is_loop = false;
	}

	stroke3DPoints.resize(0, 3);
	for(int i = 0; i < closest_vert_bindings.size()-1; i++) {
		set<int> goal;
		goal.insert(closest_vert_bindings[i + 1]);
		dijkstra_compute_paths(closest_vert_bindings[i], goal, adj_list, min_dist, previous);
		dijkstra_get_shortest_path_to(closest_vert_bindings[i + 1], previous, result_path);
	
		for(int j = result_path.size()-1; j >= 0; j--) {
			int idx = result_path[j];
			if(idx != prev) {
				added_stroke_final_vertices.push_back(idx);
				stroke3DPoints.conservativeResize(stroke3DPoints.rows() + 1, 3);
				stroke3DPoints.row(stroke3DPoints.rows() - 1) << V.row(added_stroke_final_vertices.back());
				vertex_boundary_markers[added_stroke_final_vertices.back()] = stroke_ID; //Set all vertices that are in the added stroke3DPoints to be "boundary"/constraint vertices
				prev = idx;
			}
		}
	}
	stroke3DPoints.conservativeResize(stroke3DPoints.rows() + 1, 3);
	stroke3DPoints.row(stroke3DPoints.rows() - 1) << stroke3DPoints.row(0); // Also add a copy of the first point at the end for strokes that are later added, to keep consistency with the initial stroke
	
}

void Stroke::mirror_on_backside(Eigen::VectorXi &vertex_boundary_markers, unordered_map<int, int> backside_vertex_map) {
	int original_added_stroke_size = added_stroke_final_vertices.size();
	stroke3DPoints.conservativeResize(stroke3DPoints.rows() - 1, 3); //Should remove the last point, which is a duplicate of the first, to keep consistency with original strokes

	int back_idx;
	for(int i = original_added_stroke_size - 1; i >= 0; i--) {
		back_idx = backside_vertex_map.at(added_stroke_final_vertices[i]);
		if(back_idx != added_stroke_final_vertices[i]) { //Don't add "backside" vertices for vertices that are on the original stroke
			added_stroke_final_vertices.push_back(back_idx);
			stroke3DPoints.conservativeResize(stroke3DPoints.rows() + 1, 3);
			stroke3DPoints.row(stroke3DPoints.rows() - 1) << V.row(back_idx);
			vertex_boundary_markers[added_stroke_final_vertices.back()] = stroke_ID;
			closest_vert_bindings.push_back(back_idx);
		}
	}
	stroke3DPoints.conservativeResize(stroke3DPoints.rows() + 1, 3);
	stroke3DPoints.row(stroke3DPoints.rows() - 1) << stroke3DPoints.row(0); // Also add a copy of the first point at the end for strokes that are later added, to keep consistency with the initial stroke
	is_loop = true;
}

int Stroke::get_ID() {
	return stroke_ID;
}

//Adds stroke elements at the intersection points of the original drawn stroke with mesh edges
/*void Stroke::prepare_for_cut() {
	Eigen::Matrix4f modelview = viewer.core.view * viewer.core.model;
	int faceID = -1;
	Eigen::Vector3f bc;

	int prev_p = 0;
	int start_p = prev_p;
	int next_p;
	igl::unproject_onto_mesh(stroke2DPoints.row(prev_p), modelview, viewer.core.proj, viewer.core.viewport, V, F, faceID, bc);

	int start_face = faceID;
	int n = 1;
	Eigen::RowVector3d pt(0, 0, 0);

	igl::edge_topology(V, F, EV, FE, EF);

	while(true) {
		next_p = n;
		pt = V.row(F(faceID, 0))*bc(0) + V.row(F(faceID, 1))*bc(1) + V.row(F(faceID, 2))*bc(2);

		Eigen::VectorXi forward;
		faceID = extend_path(prev_p, next_p, faceID, forward, false, modelview);
	}
}

int Stroke::extend_path(int prev_p, int next_p, int faceID, Eigen::VectorXi forward, bool front_facing_only, Eigen::Matrix4f modelview) {
	Eigen::Vector3d source, dir;
	Eigen::Vector2d tmp = stroke2DPoints.row(prev_p);
	igl::unproject_ray(tmp, modelview, viewer.core.proj, viewer.core.viewport, source, dir);
	Plane cutPlane(source, stroke3DPoints.row(prev_p), stroke3DPoints.row(next_p));

	int edge = NULL;
	pair<int, int> strokeEdge(prev_p, next_p);

	int proj_faceID = -1;
	Eigen::Vector3f bc;
	while(true) {
		if(igl::unproject_onto_mesh(stroke2DPoints.row(next_p), modelview, viewer.core.proj, viewer.core.viewport, V, F, proj_faceID, bc)){
			if(proj_faceID == faceID) {
				forward[0] = true;
				return faceID;
			}
		}

		edge = find_next_edge(strokeEdge, edge, faceID, modelview);
		if(edge == NULL) {
			return NULL;
		}

		Eigen::Vector3d v = cutPlane.cross_point(stroke3DPoints.row(EV(edge,0)), stroke3DPoints.row(EV(edge,1)));

		//TODO:: ADD PATH ELEMENT LIKE SURFACEPATH LINE 357


		faceID = (EF(edge, 0) == faceID) ? EF(edge, 1) : EF(edge, 0); //get the polygon on the other side of the edge

		if(igl::unproject_onto_mesh(stroke2DPoints.row(prev_p), modelview, viewer.core.proj, viewer.core.viewport, V, F, proj_faceID, bc)) {
			if(proj_faceID == faceID) {
				forward[0] = false;
				return faceID;
			}
		}
	}
}

int Stroke::find_next_edge(pair<int, int> strokeEdge, int prev_edge, int polygon, Eigen::Matrix4f modelview) {
	Eigen::RowVector3d start, end, stroke_start, stroke_end;
	Eigen::Vector3d tmp;
	for(int i = 0; i < 3; i++) {
		int edge = FE(polygon, i);
		if(edge != prev_edge) {
			tmp = V.row(EV(edge, 0)).transpose();
			igl::project(tmp, modelview, viewer.core.proj, viewer.core.viewport, start);

			tmp = V.row(EV(edge, 1)).transpose();
			igl::project(tmp, modelview, viewer.core.proj, viewer.core.viewport, end);

			stroke_start = stroke2DPoints.row(strokeEdge.first);
			stroke_end = stroke2DPoints.row(strokeEdge.second);
			if(edges2D_cross({stroke_start, stroke_end}, {start, end})) {
				return edge;
			}
		}
	}
	return NULL;
}


//Follows principle from https://stackoverflow.com/questions/14176776/find-out-if-2-lines-intersect but slightly different
bool Stroke::edges2D_cross(pair<Eigen::Vector2d, Eigen::Vector2d> edge1, pair<Eigen::Vector2d, Eigen::Vector2d> edge2) {
	double a0, b0, c0, a1, b1, c1;
	a0 = edge1.first[1] - edge1.second[1]; //y coordinates of start and end point of first edge
	b0 = edge1.second[0] - edge1.first[0];
	c0 = edge1.second[1] * edge1.first[0] - edge1.second[0] * edge1.first[1];
	a1 = edge2.first[1] - edge2.second[1]; //y coordinates of start and end point of first edge
	b1 = edge2.second[0] - edge2.first[0];
	c1 = edge2.second[1] * edge2.first[0] - edge2.second[0] * edge2.first[1];

	if(((a0*edge2.first[0] + b0*edge2.first[1] + c0)*(a0*edge2.second[0] + b0*edge2.second[1] + c0) <= 0) &&
		((a1*edge1.first[0] + b1*edge1.first[1] + c1)*(a1*edge1.second[0] + b1*edge1.second[1] + c1) <= 0)) {
		return true;
	} else {
		return false;
	}
}
*/

Eigen::MatrixXd Stroke::get_V() const {
	return V;
}

Eigen::MatrixXi Stroke::get_F() const{
	return F;
}

Eigen::MatrixX2d Stroke::get_stroke2DPoints() const {
	return stroke2DPoints;
}