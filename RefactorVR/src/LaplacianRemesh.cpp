#include "LaplacianRemesh.h"
#include "CleanStroke3D.h"
#include <igl/edge_topology.h>
#include <igl/vertex_triangle_adjacency.h>
#include <igl/slice.h>
#include <igl/adjacency_list.h>
#include <igl/project.h>
#include <igl/cat.h>


using namespace std;
using namespace igl;

bool LaplacianRemesh::is_front_loop = true; //Determines extrusion or cut
bool LaplacianRemesh::remove_inside_faces = true;

Eigen::MatrixXi LaplacianRemesh::EV, LaplacianRemesh::FE, LaplacianRemesh::EF;
vector<vector<int>> LaplacianRemesh::VV;

Eigen::VectorXi LaplacianRemesh::remesh_cut_remove_inside(Mesh & m, SurfacePath & surface_path, Eigen::Matrix4f model, Eigen::Matrix4f view, Eigen::Matrix4f proj, Eigen::Vector4f viewport, bool& remesh_success, int cut_clicked_face) {
	is_front_loop = false;
	remove_inside_faces = true;
	adjacency_list(m.F, VV);
	return remesh(m, surface_path, model, view, proj, viewport, remesh_success, cut_clicked_face);
}

Eigen::VectorXi LaplacianRemesh::remesh_extrusion_remove_inside(Mesh & m, SurfacePath & surface_path, Eigen::Matrix4f model, Eigen::Matrix4f view, Eigen::Matrix4f proj, Eigen::Vector4f viewport, bool& remesh_success) {
	is_front_loop = true;
	remove_inside_faces = true;
	adjacency_list(m.F, VV);
	return remesh(m, surface_path, model, view, proj, viewport, remesh_success);
}

bool LaplacianRemesh::remesh_open_path(Mesh& m, Stroke& open_path_stroke) {
	SurfacePath surface_path;
	surface_path.create_from_open_path(open_path_stroke);

	Eigen::MatrixXi start_F = m.F;
	Eigen::MatrixXd start_V = m.V;
	Eigen::VectorXi start_part_of_original_stroke = m.part_of_original_stroke;
	Eigen::VectorXi start_vertex_boundary_markers = m.vertex_boundary_markers;
	Eigen::VectorXi start_sharp_edge = m.sharp_edge;
	Eigen::VectorXi start_new_mapped_indices = m.new_mapped_indices;

	vector<bool> dirty_face(m.F.rows());
	Eigen::VectorXi dirty_vertices(m.V.rows());
	dirty_vertices.setZero();

	int face = -1, edge = -1;
	vector<PathElement> path = surface_path.get_path();

	m.new_mapped_indices.resize(m.V.rows()); //Resize the map from old to new (clean) vertex indices to allow it to contain the number of vertices that are in the mesh at the start

	vector<int> sharp_edge_indices;
	for (int i = 0; i < m.sharp_edge.rows(); i++) {
		if (m.sharp_edge[i]) {
			sharp_edge_indices.push_back(i);
		}
	}

	igl::edge_topology(m.V, m.F, EV, FE, EF);
	adjacency_list(m.F, VV);

	Eigen::MatrixXi sharpEV;
	Eigen::VectorXi sharpEV_row_idx, sharpEV_col_idx(2);
	sharpEV_row_idx = Eigen::VectorXi::Map(sharp_edge_indices.data(), sharp_edge_indices.size());
	sharpEV_col_idx.col(0) << 0, 1;
	igl::slice(EV, sharpEV_row_idx, sharpEV_col_idx, sharpEV); //Keep only the sharp edges in the original mesh

	for (int i = 0; i < path.size(); i++) {
		if (path[i].get_type() == PathElement::FACE) {
			face = path[i].get_ID();
			dirty_vertices[m.F(face, 0)] = true; //These will ensure that all vertices of the end faces will also be in outer_boundary_vertices
			dirty_vertices[m.F(face, 1)] = true;
			dirty_vertices[m.F(face, 2)] = true;

		}
		else {
			edge = path[i].get_ID();
			dirty_vertices[EV(edge, 0)] += 1;
			dirty_vertices[EV(edge, 1)] += 1;
			dirty_face[EF(edge, 0)] = true;
			dirty_face[EF(edge, 1)] = true;
			face = (EF(edge, 0) == face) ? EF(edge, 1) : EF(edge, 0); //get the polygon on the other side of the edge
		}
	}

	//Collect vertices on the boundary
	vector<int> outer_boundary_vertices;
	for (int i = 0; i < m.V.rows(); i++) {
		if (dirty_vertices[i] > 0) {
			outer_boundary_vertices.push_back(i);
		}
		m.new_mapped_indices[i] = i;
	}

	//Remove dirty faces
	vector<int> clean_faces;
	for (int i = 0; i < m.F.rows(); i++) {
		if (!dirty_face[i]) {
			clean_faces.push_back(i);
		}
	}

	Eigen::MatrixXi tmp_F;
	Eigen::MatrixXi prev_F = m.F;
	Eigen::VectorXi row_idx, col_idx(3);
	row_idx = Eigen::VectorXi::Map(clean_faces.data(), clean_faces.size());
	col_idx.col(0) << 0, 1, 2;
	igl::slice(m.F, row_idx, col_idx, tmp_F); //Keep only the clean faces in the mesh
	m.F = tmp_F;

	//NOTE: Output from sort_boundary_vertices is not necessarily counter-clockwise
	//TODO: Handle errors nicely here (e.g. make this method return a bool that gets checked)
	try {
		for (int i = 0; i < outer_boundary_vertices.size(); i++) {
			std::cout << outer_boundary_vertices[i] << "   ";
		}
		std::cout << std::endl;
		outer_boundary_vertices = sort_boundary_vertices(path[0].get_vertex(), outer_boundary_vertices, m);
	}
	catch (int ex) {
		if (ex == -1) {
			std::cerr << "Could not process adding stroke. Please try again. " << std::endl;
			m.F = prev_F;
			//remesh_success = false;
			return false;
		}
		else {
			std::cout << "Something weird happened. Check what's going on." << std::endl;
			m.F = prev_F;
			//remesh_success = false;
			return false;
		}
	}

	Eigen::MatrixXi replacing_sharp_edges(0, 2);
	for (int i = 0; i < path.size(); i++) {
		if (path[i].get_type() == PathElement::EDGE && m.sharp_edge[path[i].get_ID()]) {
			//TODO: store that we are removing/splitting a sharp edge, so that we can later set the 2 resulting new edges to be sharp
			path[i].fixed = true;
			replacing_sharp_edges.conservativeResize(replacing_sharp_edges.rows() + 2, Eigen::NoChange);
			replacing_sharp_edges.row(replacing_sharp_edges.rows() - 2) << m.V.rows() + i, EV(path[i].get_ID(), 0); //Pair of original edge start index and middle vertex index
			replacing_sharp_edges.row(replacing_sharp_edges.rows() - 1) << m.V.rows() + i, EV(path[i].get_ID(), 1); //Pair of original edge end index and middle vertex index
		}
		else {
			path[i].fixed = false;
		}
	}



	Eigen::RowVector3d mean_viewpoint = compute_mean_viewpoint(m, outer_boundary_vertices);
	Eigen::Matrix4f modelview = open_path_stroke.viewer.oculusVR.get_start_action_view() * open_path_stroke.viewer.core.get_model();

	double unit_length = compute_average_distance_between_onPolygon_vertices(path);
	Eigen::MatrixXd resampled_path = CleanStroke3D::resample_by_length_with_fixes(path, unit_length);
	if (!is_counter_clockwise_boundaries(resampled_path, modelview, open_path_stroke.viewer.core.get_proj(), open_path_stroke.viewer.core.viewport, mean_viewpoint, true)) {
		resampled_path = resampled_path.colwise().reverse().eval();
	}

	vector<int> path_vertices;
	int original_V_size = m.V.rows(); //TODO: Maybe have to redo below that it doesn't take the last vertex from the resampled path (probably not since it's not a loop)
	for (int i = 0; i < resampled_path.rows(); i++) {
		path_vertices.push_back(original_V_size + i); //Store the vertex indices (in m.V) of the path vertices
	}

	vector<PathElement> new_surface_path;
	for (int i = 0; i < resampled_path.rows(); i++) { //Make the new surfacePath looped (unlike the path_vertices indices) //TODO: remove this comment? (open path ain't looped)
		new_surface_path.push_back(PathElement(resampled_path.row(i), m.V.rows() + i));
	}
	surface_path.set_path(new_surface_path);

	update_mesh_values(m, resampled_path, surface_path.get_origin_stroke_ID(), m.V.rows(), false);




	//Add interior path vertices in reversed direction to the path, so we can complete the loop around it
	for (int i = path_vertices.size() - 2; i > 0; i--) {
		path_vertices.push_back(path_vertices[i]);
	}

	//Make the outer_boundary_vertices CCW
	Eigen::MatrixXd tmp_V;
	Eigen::VectorXi row_idx2, col_idx2(3);
	row_idx2 = Eigen::VectorXi::Map(outer_boundary_vertices.data(), outer_boundary_vertices.size());
	col_idx2.col(0) << 0, 1, 2;
	igl::slice(m.V, row_idx2, col_idx2, tmp_V); //Keep only the clean "boundary" vertices in the mesh

	if (!is_counter_clockwise_boundaries(tmp_V, modelview, open_path_stroke.viewer.core.get_proj(), open_path_stroke.viewer.core.viewport, mean_viewpoint, !is_front_loop)) {
		reverse(outer_boundary_vertices.begin(), outer_boundary_vertices.end());
	}

	//TODO: for open paths we may need to add the extra elseif statement to stitch(), but not sure if this will work together with stitching for cut (since this is also looped over to a backside)
	stitch(path_vertices, outer_boundary_vertices, m); //TODO: need to fix that stitching sometimes fails


	Eigen::MatrixXi sharpEVcat = igl::cat(1, sharpEV, replacing_sharp_edges);
	try {
		update_sharp_edges(m, sharpEVcat);
	}
	catch (int ex) {
		if (ex == -1) {
			std::cerr << "Remeshing results in a non edge-manifold mesh, which cannot be processed. Please try again. " << std::endl;
			m.F = start_F;
			m.V = start_V;
			m.new_mapped_indices = start_new_mapped_indices;
			m.part_of_original_stroke = start_part_of_original_stroke;
			m.vertex_boundary_markers = start_vertex_boundary_markers;
			m.sharp_edge = start_sharp_edge;
			return false;
		}
	}


	//TODO: Maybe add  MeshExtrusion::post_extrude_main_update_points(Stroke &stroke, Eigen::MatrixXd new_positions) here

	update_open_path_points_and_bindings(open_path_stroke, surface_path);
	return true;
}

bool LaplacianRemesh::remesh_cutting_path(Mesh& m, Stroke& cut_path_stroke) {
	remove_inside_faces = false;

	return true;
}

Eigen::VectorXi LaplacianRemesh::remesh(Mesh& m, SurfacePath& surface_path, Eigen::Matrix4f model, Eigen::Matrix4f view, Eigen::Matrix4f proj, Eigen::Vector4f viewport, bool& remesh_success, int cut_clicked_face) {
	Eigen::MatrixXi start_F = m.F;
	Eigen::MatrixXd start_V = m.V;
	Eigen::VectorXi start_part_of_original_stroke = m.part_of_original_stroke;
	Eigen::VectorXi start_vertex_boundary_markers = m.vertex_boundary_markers;
	Eigen::VectorXi start_sharp_edge = m.sharp_edge;
	Eigen::VectorXi start_new_mapped_indices = m.new_mapped_indices;

	vector<bool> dirty_face(m.F.rows());
	Eigen::VectorXi dirty_vertices(m.V.rows());
	dirty_vertices.setZero();
	m.new_mapped_indices.resize(m.V.rows()); //Resize the map from old to new (clean) vertex indices to allow it to contain the number of vertices that are in the mesh at the start

	vector<int> sharp_edge_indices;
	for (int i = 0; i < m.sharp_edge.rows(); i++) {
		if (m.sharp_edge[i]) {
			sharp_edge_indices.push_back(i);
		}
	}

	vector<vector<int>> VF, VI;
	igl::edge_topology(m.V, m.F, EV, FE, EF);
	igl::vertex_triangle_adjacency(m.V.rows(), m.F, VF, VI);
	Eigen::MatrixXi startEV = EV;
	Eigen::MatrixXd startV = m.V;

	Eigen::MatrixXi sharpEV;
	Eigen::VectorXi sharpEV_row_idx, sharpEV_col_idx(2);
	sharpEV_row_idx = Eigen::VectorXi::Map(sharp_edge_indices.data(), sharp_edge_indices.size());
	sharpEV_col_idx.col(0) << 0, 1;
	igl::slice(EV, sharpEV_row_idx, sharpEV_col_idx, sharpEV); //Keep only the sharp edges in the original mesh

	int face = -1, edge = -1;
	vector<PathElement> path = surface_path.get_path();
	if (!is_front_loop) { //Append first element at the end in case of cutting, to make it a loop and compatible with extrusion's loop
		path.push_back(path[0]);
	}

	//Find outside and inside vertices. +1 means that the vertex will stay, -1 means it will be destroyed
	for (int i = 0; i < path.size(); i++) {
		if (path[i].get_type() == PathElement::FACE) {
			face = path[i].get_ID();
		}
		else {
			edge = path[i].get_ID();
			if (EF(edge, 0) == face) {
				dirty_vertices[EV(edge, 0)] += 1;
				dirty_vertices[EV(edge, 1)] -= 1;
			}
			else {
				dirty_vertices[EV(edge, 0)] -= 1;
				dirty_vertices[EV(edge, 1)] += 1;
			}
			face = (EF(edge, 0) == face) ? EF(edge, 1) : EF(edge, 0); //get the polygon on the other side of the edge
		}
	}

	//Collect vertices on the boundary
	vector<int> outer_boundary_vertices;
	vector<int> inner_boundary_vertices;
	for (int i = 0; i < m.V.rows(); i++) {
		if (dirty_vertices[i] < 0) { //vertices that are removed
			inner_boundary_vertices.push_back(i);
		}
		else if (dirty_vertices[i] > 0) { //vertices that stay
			outer_boundary_vertices.push_back(i);
		}
	}

	if (outer_boundary_vertices.size() == 0 || inner_boundary_vertices.size() == 0) { //There are no vertices left behind (cut that removes everthing or extrude that doesn't include at least 1 vertex) or none are removed
		std::cerr << "This action either removes all or no vertices, which is not possible. Please try again. " << std::endl;
		remesh_success = false;
		return Eigen::VectorXi::Zero(1);
	}

	//Collect faces along the path
	for (int i = 0; i < path.size(); i++) {
		if (path[i].get_type() == PathElement::EDGE) {
			dirty_face[EF(path[i].get_ID(), 0)] = true;
			dirty_face[EF(path[i].get_ID(), 1)] = true;
		}
	}

	if (remove_inside_faces) {
		for (int i = 0; i < m.V.rows(); i++) {
			if (dirty_vertices[i] < 0) {
				for (int j = 0; j < VF[i].size(); j++) {
					if (!dirty_face[VF[i][j]]) {
						propagate_dirty_faces(VF[i][j], dirty_face);
					}
				}
			}
		}
	}


	//Revert dirty/clean faces and vertices in case the user clicked on the "clean" part when cutting
	if (cut_clicked_face >= 0 && !dirty_face[cut_clicked_face] && remove_inside_faces) {
		for (int i = 0; i < dirty_face.size(); i++) {
			dirty_face[i] = !dirty_face[i];
		}

		dirty_vertices *= -1;

		inner_boundary_vertices.swap(outer_boundary_vertices);

		//Collect faces along the path
		for (int i = 0; i < path.size(); i++) {
			if (path[i].get_type() == PathElement::EDGE) {
				dirty_face[EF(path[i].get_ID(), 0)] = true;
				dirty_face[EF(path[i].get_ID(), 1)] = true;
			}
		}
	}

	//Remove dirty faces
	vector<int> clean_faces;
	for (int i = 0; i < m.F.rows(); i++) {
		if (!dirty_face[i]) {
			clean_faces.push_back(i);
		}
	}


	if (clean_faces.size() == 0) { //There are no faces left behind (cut that removes everything or extrude that doesn't include at least 1 vertex)
		if (cut_clicked_face == -1) { //Extrusion
			std::cerr << "The base does not include at least 1 vertex.  Please try again." << std::endl;
		}
		else { //Cut
			std::cerr << "This removes all vertices. Please try again." << std::endl;
		}
		remesh_success = false;
		return Eigen::VectorXi::Zero(1);
	}

	//Determine which vertices are clean
	Eigen::VectorXi vertex_is_clean = Eigen::VectorXi::Zero(m.V.rows());
	vector<int> clean_vertices;
	for (int i = 0; i < m.V.rows(); i++) {
		if (dirty_vertices[i] == 0) { //Vertices with a - value are already known to be removed and vertices with a + value are already known to stay
			for (int j = 0; j < VF[i].size(); j++) {
				if (!dirty_face[VF[i][j]]) { //If the vertex is part of a clean face, then the vertex has to be clean
					vertex_is_clean(i) = 1;
					clean_vertices.push_back(i);
					break; //Can skip checking all other adjacent faces for this vertex
				}
			}
		}
		else if (dirty_vertices[i] > 0) {
			vertex_is_clean(i) = 1;
			clean_vertices.push_back(i);
		}
	}

	Eigen::MatrixXi tmp_F;
	Eigen::VectorXi row_idx, col_idx(3);
	row_idx = Eigen::VectorXi::Map(clean_faces.data(), clean_faces.size());
	col_idx.col(0) << 0, 1, 2;
	igl::slice(m.F, row_idx, col_idx, tmp_F); //Keep only the clean faces in the mesh
	m.F = tmp_F;

	//NOTE: Output from sort_boundary_vertices is not necessarily counter-clockwise
	try {
		outer_boundary_vertices = sort_boundary_vertices(path[0].get_vertex(), outer_boundary_vertices, m);
	}
	catch (int ex) {
		if (ex == -1) {
			std::cerr << "Could not process cut stroke or extrusion base. Please try again. " << std::endl;
			m.F = start_F;
			remesh_success = false;
			return Eigen::VectorXi::Zero(1);
		}
		else {
			std::cout << "Something weird happened. Check what's going on." << std::endl;
			m.F = start_F;
			remesh_success = false;
			return Eigen::VectorXi::Zero(1);
		}
	}

	Eigen::RowVector3d mean_viewpoint = compute_mean_viewpoint(m, inner_boundary_vertices);

	Eigen::MatrixXd tmp_V;
	Eigen::VectorXi tmp_part_of, tmp_markers;
	int size_before_removing = m.V.rows();
	row_idx = Eigen::VectorXi::Map(clean_vertices.data(), clean_vertices.size());
	col_idx.col(0) << 0, 1, 2;
	igl::slice(m.V, row_idx, col_idx, tmp_V);
	m.V = tmp_V; //Takes care of removing vertices that no longer exist due to all their adjacent faces being removed

	col_idx.resize(1);
	col_idx.col(0) << 0;
	igl::slice(m.part_of_original_stroke, row_idx, col_idx, tmp_part_of);
	igl::slice(m.vertex_boundary_markers, row_idx, col_idx, tmp_markers);
	m.part_of_original_stroke = tmp_part_of; //Compact part_of_original_stroke by removing the values for vertices that are now removed
	m.vertex_boundary_markers = tmp_markers; //Compact vertex_boundary_markers by removing the values for vertices that are now removed


	//Map indices in m.V at the start of remeshing to indices in m.V after remeshing
	int count = 0;
	m.new_mapped_indices.setConstant(m.new_mapped_indices.rows(), -1); //Initialize all values to -1
	for (int i = 0; i < vertex_is_clean.rows(); i++) {
		if (vertex_is_clean[i]) {
			m.new_mapped_indices[i] = count;
			count++;
		}
	}

	update_face_indices(m);

	//Up to here disappearing vertices and faces have been removed from the mesh and corresponding indices/markers have been updated
	Eigen::MatrixXi replacing_sharp_edges(0, 2);
	for (int i = 0; i < path.size(); i++) {
		if (path[i].get_type() == PathElement::EDGE && m.sharp_edge[path[i].get_ID()]) {
			//TODO: store that we are removing/splitting a sharp edge, so that we can later set the 2 resulting new edges to be sharp
			path[i].fixed = true;
			replacing_sharp_edges.conservativeResize(replacing_sharp_edges.rows() + 2, Eigen::NoChange);
			replacing_sharp_edges.row(replacing_sharp_edges.rows() - 2) << vertex_is_clean.rows() + i, startEV(path[i].get_ID(), 0); //Pair of original edge start index and middle vertex index
			replacing_sharp_edges.row(replacing_sharp_edges.rows() - 1) << vertex_is_clean.rows() + i, startEV(path[i].get_ID(), 1); //Pair of original edge end index and middle vertex index
		}
		else {
			path[i].fixed = false;
		}
	}

	double unit_length = (is_front_loop) ? compute_average_distance_between_onPolygon_vertices(path) : compute_average_length_of_crossing_edges(path, startV, startEV);
	if (CleanStroke3D::get_stroke_length(path, 0, path.size() - 1) / unit_length < 12) { //We'd end up with less than 12 samples
		unit_length = CleanStroke3D::get_stroke_length(path, 0, path.size() - 1) / 12; //Adapt length such that we get 12 samples
	}

	Eigen::MatrixXd resampled_path = CleanStroke3D::resample_by_length_with_fixes(path, unit_length);
	Eigen::Matrix4f modelview = view * model;
	if (!is_counter_clockwise_boundaries(resampled_path, modelview, proj, viewport, mean_viewpoint, true)) {
		resampled_path = resampled_path.colwise().reverse().eval();
	}


	vector<int> path_vertices;
	int original_V_size = m.V.rows();
	for (int i = 0; i < resampled_path.rows() - 1; i++) { //Do not use the last path vertex, as it is a copy of the first and creates unwanted behaviour
		path_vertices.push_back(original_V_size + i); //Store the vertex indices (in m.V) of the path vertices
	}

	update_mesh_values(m, resampled_path, surface_path.get_origin_stroke_ID(), vertex_is_clean.rows(), true);

	vector<PathElement> new_surface_path;
	for (int i = 0; i < resampled_path.rows(); i++) { //Make the new surfacePath looped (unlike the path_vertices indices)
		new_surface_path.push_back(PathElement(resampled_path.row(i), size_before_removing + i));
	}
	surface_path.set_path(new_surface_path);

	//TODO: need to make sure that we don't add sharp edges for added strokes (that are not a cut or extrusion)
	Eigen::MatrixXi added_sharpEV(resampled_path.rows() - 1, 2);
	for (int i = 0; i < resampled_path.rows() - 1; i++) {
		added_sharpEV.row(i) << vertex_is_clean.rows() + i, vertex_is_clean.rows() + ((i + 1) + resampled_path.rows() - 1) % (resampled_path.rows() - 1); //Use mod (path.size() - 1) to compensate for start point that is in path twice
	}


	//Update the outer_boundary_vertices (which are indices into V before it was sliced to keep only the clean vertices)
	for (int i = 0; i < outer_boundary_vertices.size(); i++) {
		outer_boundary_vertices[i] = m.new_mapped_indices[outer_boundary_vertices[i]];
	}

	//Make the outer_boundary_vertices CCW
	Eigen::VectorXi row_idx2, col_idx2(3);
	row_idx2 = Eigen::VectorXi::Map(outer_boundary_vertices.data(), outer_boundary_vertices.size());
	col_idx2.col(0) << 0, 1, 2;
	igl::slice(m.V, row_idx2, col_idx2, tmp_V); //Keep only the clean "boundary" vertices in the mesh

	if (!is_counter_clockwise_boundaries(tmp_V, modelview, proj, viewport, mean_viewpoint, !is_front_loop)) {
		reverse(outer_boundary_vertices.begin(), outer_boundary_vertices.end());
		if (!remove_inside_faces) {
			reverse(inner_boundary_vertices.begin(), inner_boundary_vertices.end());
		}
	}

	stitch(path_vertices, outer_boundary_vertices, m); //TODO: need to fix that stitching sometimes fails
	if (!remove_inside_faces) {
		for (int i = 0; i < inner_boundary_vertices.size(); i++) {
			inner_boundary_vertices[i] = m.new_mapped_indices[inner_boundary_vertices[i]];
		}
		stitch(path_vertices, inner_boundary_vertices, m);
	}

	Eigen::MatrixXi sharpEVcat;
	if (remove_inside_faces) { //True for extrusion and cut, who both create a sharp curve along their path
		sharpEVcat = igl::cat(1, sharpEV, added_sharpEV); //Add sharp edges that are created due to the newly added curve
		sharpEVcat = igl::cat(1, sharpEVcat, replacing_sharp_edges); //Add sharp edges that are created due to old sharp edges being broken
	}
	else {
		sharpEVcat = igl::cat(1, sharpEV, replacing_sharp_edges); //Add sharp edges that are created due to old sharp edges being broken
	}


	try {
		update_sharp_edges(m, sharpEVcat);
	}
	catch (int ex) {
		if (ex == -1) {
			std::cerr << "Remeshing results in a non edge-manifold mesh, which cannot be processed. Please try again. " << std::endl;
			m.F = start_F;
			m.V = start_V;
			m.new_mapped_indices = start_new_mapped_indices;
			m.part_of_original_stroke = start_part_of_original_stroke;
			m.vertex_boundary_markers = start_vertex_boundary_markers;
			m.sharp_edge = start_sharp_edge;
			remesh_success = false;
			return Eigen::VectorXi::Zero(1);
		}
	}
	return Eigen::VectorXi::Map(path_vertices.data(), path_vertices.size());
}

/** Compute the mean of the boundary vertices on the side of the mesh that is removed.
	In the case of CUT this is used as a viewpoint to determine the orientation (CW or CCW) of the boundary vertices that remain.
	Do this before slicing m.V and removing the vertices that we need. **/
Eigen::RowVector3d LaplacianRemesh::compute_mean_viewpoint(Mesh &m, vector<int> inner_boundary_vertices) {
	Eigen::MatrixXd removed_V;
	Eigen::VectorXi row_idx_removed = Eigen::VectorXi::Map(inner_boundary_vertices.data(), inner_boundary_vertices.size());
	Eigen::VectorXi col_idx_removed(3);
	col_idx_removed.col(0) << 0, 1, 2;
	igl::slice(m.V, row_idx_removed, col_idx_removed, removed_V);

	return removed_V.colwise().mean();
}

/** Updates the faces with new vertex indexing. **/
void LaplacianRemesh::update_face_indices(Mesh &m) {
	for (int i = 0; i < m.F.rows(); i++) {
		for (int j = 0; j < 3; j++) {
			m.F(i, j) = m.new_mapped_indices[m.F(i, j)];
		}
	}
}

/** Adds values for a bunch of tracking variables for the new mesh vertices. **/
void LaplacianRemesh::update_mesh_values(Mesh& m, Eigen::MatrixXd path, int stroke_ID, int new_mapped_start, bool hold_back_due_to_loop) {
	int size_before = m.V.rows();
	m.V.conservativeResize(m.V.rows() + path.rows() - hold_back_due_to_loop, Eigen::NoChange);
	m.part_of_original_stroke.conservativeResize(m.part_of_original_stroke.rows() + path.rows() - hold_back_due_to_loop);
	m.vertex_boundary_markers.conservativeResize(m.vertex_boundary_markers.rows() + path.rows() - hold_back_due_to_loop);
	m.new_mapped_indices.conservativeResize(m.new_mapped_indices.rows() + path.rows() - hold_back_due_to_loop, Eigen::NoChange);
	for (int i = 0; i < path.rows() - hold_back_due_to_loop; i++) {
		m.V.row(size_before + i) << path.row(i);
		m.part_of_original_stroke[size_before + i] = 0;
		m.vertex_boundary_markers[size_before + i] = stroke_ID;
		m.new_mapped_indices(new_mapped_start + i) = size_before + i;
	}
}

/** Updates the old sharp edges after the mesh topology changed and inserts newly generated sharp edges. **/
void LaplacianRemesh::update_sharp_edges(Mesh& m, Eigen::MatrixXi& all_sharpEV) {
	if (!igl::is_edge_manifold(m.F)) {
		throw -1;
		return;
	}
	igl::edge_topology(m.V, m.F, EV, FE, EF);
	m.sharp_edge.resize(EV.rows());
	m.sharp_edge.setZero();

	int start, end, equal_pos;
	Eigen::VectorXi col1Equals, col2Equals;
	for (int i = 0; i < all_sharpEV.rows(); i++) {
		std::cout << all_sharpEV(i, 0) << "   " << all_sharpEV(i, 1) << std::endl;
		start = m.new_mapped_indices(all_sharpEV(i, 0));
		end = m.new_mapped_indices(all_sharpEV(i, 1));
		if (start == -1 || end == -1) { //Sharp edge no longer exists
			continue;
		}

		col1Equals = EV.col(0).cwiseEqual(min(start, end)).cast<int>();
		col2Equals = EV.col(1).cwiseEqual(max(start, end)).cast<int>();
		(col1Equals + col2Equals).maxCoeff(&equal_pos); //Find the row that contains both vertices of this edge

		m.sharp_edge[equal_pos] = 1; //Set this edge to be sharp
	}
}

/** Floodfills all dirty faces. **/
void LaplacianRemesh::propagate_dirty_faces(int face, vector<bool>& dirty_face) {
	dirty_face[face] = true;
	for (int i = 0; i < 3; i++) {
		int next_face = (EF(FE(face, i), 0) == face) ? EF(FE(face, i), 1) : EF(FE(face, i), 0); //get the adjacent polygon
		if (!dirty_face[next_face]) {
			propagate_dirty_faces(next_face, dirty_face);
		}
	}
}

//TODO: possible way to improve this: iterate over the surface_path, and extract the edges that it crosses (which are in order). For each edge look at it's adjacent vertices and keep the one that is on the "staying" side. These vertices should now also be in order.
//Then rotate this ordered sequence around so it starts with the vertex closest to start_vertex.
/** Sorts the boundary vertices such that they make up a continuous path that starts with the vertex that's closest to start_vertex.
	Takes care of backtracking over original stroke boundary vertices that is sometimes necessary. **/
vector<int> LaplacianRemesh::sort_boundary_vertices(Eigen::Vector3d start_vertex, vector<int> boundary_vertices, Mesh& m) {
	igl::edge_topology(m.V, m.F, EV, FE, EF);

	vector<int> sorted_boundary_vertices;
	int start_v = find_closest(boundary_vertices, start_vertex, m);
	int v = start_v;
	int prev = start_v;
	int iter = 0;

	Eigen::VectorXi tmp1, tmp2;
	int min_idx, max_idx, equal_pos, v_pos, max_val, backtracking_option_idx, disconnected_option_idx;
	bool backtracking_option_available = false, disconnected_option_available = false;
	while (true) {
		iter++;
		if (iter == 10 * boundary_vertices.size()) { //We're not finding a solution
			throw - 1;
			return sorted_boundary_vertices; //Will be discarded anyway
		}

		sorted_boundary_vertices.push_back(v);
		for (int i = 0; i < VV[v].size(); i++) {
			if (std::find(boundary_vertices.begin(), boundary_vertices.end(), VV[v][i]) != boundary_vertices.end()) { //Only continue if the adjacent vertex is also a boundary vertex
				if (v < VV[v][i]) {
					min_idx = v;
					v_pos = 0;
					max_idx = VV[v][i];
				}
				else {
					min_idx = VV[v][i];
					v_pos = 1;
					max_idx = v;
				}
				tmp1 = EV.col(0).cwiseEqual(min_idx).cast<int>();
				tmp2 = EV.col(1).cwiseEqual(max_idx).cast<int>();
				max_val = (tmp1 + tmp2).maxCoeff(&equal_pos); //Tmp1 and tmp2 will contain a 1 on the positions where they equal the min and max_idx. When adding them, we'll get 2 on the row that contains both of them
				if (max_val == 2 && EF(equal_pos, v_pos) == -1) { //An edge between these vertices still exists and has a NULL face on the correct side of the edge, so proceed to adjacent vertex
					prev = v;
					v = VV[v][i];
					disconnected_option_available = false; //Remove the disconnected edge option
					backtracking_option_available = false; //Remove the backtracking option
					break;
				}
				else if (max_val != 2 && VV[v][i] != prev) { //At least one of the vertices was only connected to faces that have been removed (e.g. the edge was along the mesh boundary and is cut off)
					//First check if we can find an adjacent vertex for which a connecting edge still exists, otherwise choose the current vertex
					disconnected_option_idx = VV[v][i]; //Remember the vertex along a disconnected edge that we could fall back to
					disconnected_option_available = true;
					continue;
				}
				else if (max_val != 2 && VV[v][i] == prev && !disconnected_option_available) { //Backtracking over a disconnected edge (only set this option when we're out of non-backtracking disconnected options
					backtracking_option_idx = VV[v][i]; //Remember the vertex that we could backtrack to
					backtracking_option_available = true;
					continue;
				}
			}
		}

		if (disconnected_option_available) { //If we still have an adjacent vertex along a disconnected edge buffered, use it.
			prev = v;
			v = disconnected_option_idx;
			disconnected_option_available = false;
		}
		else if (backtracking_option_available) { //If we still have the backtracking vertex buffered, use it.
			prev = v;
			v = backtracking_option_idx;
			backtracking_option_available = false;
		}

		if (v == start_v) {
			break;
		}
	}
	return sorted_boundary_vertices;
}

/** Finds the closest vertex to point base. **/
int LaplacianRemesh::find_closest(vector<int> vertices, Eigen::Vector3d base, Mesh& m) {
	int closest = vertices[0];
	double min = (m.V.row(vertices[0]) - base.transpose()).norm();
	double d;

	for (int i = 1; i < vertices.size(); i++) {
		d = (m.V.row(vertices[i]) - base.transpose()).norm();
		if (d < min) {
			min = d;
			closest = vertices[i];
		}
	}
	return closest;
}

/** Stitches the boundary vertices (that are old, already existing mesh vertices) to the path vertices of the newly drawn stroke. **/
void LaplacianRemesh::stitch(std::vector<int> path_vertices, std::vector<int> boundary_vertices, Mesh& m) {
	int path_idx = 0, outer_idx = 0;
	boundary_vertices = reorder(boundary_vertices, m.V.row(path_vertices[0]), m);

	Eigen::RowVector3d start_path_v = m.V.row(path_vertices[0]), start_outer_v = m.V.row(boundary_vertices[0]);
	Eigen::RowVector3d path_v = start_path_v, outer_v = start_outer_v, next_path_v, next_outer_v;
	int path_v_idx = path_vertices[0], outer_v_idx = boundary_vertices[0], next_path_v_idx, next_outer_v_idx;

	bool proceed_outer_v, prev_proceed_path = false;

	while (true) {
		next_path_v = m.V.row(path_vertices[(path_idx + 1) % path_vertices.size()]);
		next_outer_v = m.V.row(boundary_vertices[(outer_idx + 1) % boundary_vertices.size()]);
		next_path_v_idx = path_vertices[(path_idx + 1) % path_vertices.size()];
		next_outer_v_idx = boundary_vertices[(outer_idx + 1) % boundary_vertices.size()];

		proceed_outer_v = true;
		if (outer_idx == boundary_vertices.size()) {
			proceed_outer_v = false;
		}
		else if (path_idx == path_vertices.size()) {
			proceed_outer_v = true;
		}
		else if (path_idx == path_vertices.size() / 2 && prev_proceed_path == true) { //The first thing after we created a triangle using the "last" path vertex ALWAYS has to be moving the outer 
			std::cout << "Do this. " << std::endl;
			proceed_outer_v = true;
		}
		else if ((path_v - next_path_v).dot(outer_v - next_outer_v) < 0 && path_vertices.size() / 2 == path_idx) {		//TODO: add extra else if for "open paths" (added curves?)
			std::cout << path_idx << "     " << path_vertices.size() / 2 << std::endl;
			//	if (path_vertices.size() / 2 == path_idx) {
			proceed_outer_v = true;
			std::cout << " because of idx" << std::endl;

			/*	}
				else if ((path_v - next_outer_v).norm() < (outer_v - next_path_v).norm()) {
					proceed_outer_v = true;
					std::cout << " because of dist" << std::endl;
				}
				else {
					proceed_outer_v = false;
				}*/
		}
		else if ((path_v - next_outer_v).norm() < (outer_v - next_path_v).norm()) {
			proceed_outer_v = true;
		}
		else {
			proceed_outer_v = false;
		}

		//Add faces
		//TODO: maybe need to add teddy.LapLacianRemesh' special handling of when we're at the start again
		if (proceed_outer_v) {
			m.F.conservativeResize(m.F.rows() + 1, Eigen::NoChange);
			m.F.row(m.F.rows() - 1) << path_v_idx, outer_v_idx, next_outer_v_idx;
			std::cout << "Added outer: " << path_v_idx << "  " << outer_v_idx << "   " << next_outer_v_idx << std::endl;
			std::cout << "Coor: " << path_v << std::endl << outer_v << std::endl << next_outer_v << std::endl << std::endl;
			outer_v = next_outer_v;
			outer_v_idx = next_outer_v_idx;
			outer_idx++;
			prev_proceed_path = false;
		}
		else {
			m.F.conservativeResize(m.F.rows() + 1, Eigen::NoChange);
			m.F.row(m.F.rows() - 1) << next_path_v_idx, path_v_idx, outer_v_idx;
			std::cout << "Added path: " << next_path_v_idx << "  " << path_v_idx << "   " << outer_v_idx << std::endl;
			std::cout << "Coor: " << next_path_v << std::endl << path_v << std::endl << outer_v << std::endl << std::endl;
			path_v = next_path_v;
			path_v_idx = next_path_v_idx;
			path_idx++;
			prev_proceed_path = true;
		}

		if (path_idx == path_vertices.size() && outer_idx == boundary_vertices.size()) {
			break;
		}
	}
}

/** Reorders the boundary vertices so they start from the vertex closest to start_v. Keeps their original order **/
vector<int> LaplacianRemesh::reorder(vector<int> boundary_vertices, Eigen::Vector3d start_v, Mesh& m) {
	int index = 0;
	double min = INFINITY;
	for (int i = 0; i < boundary_vertices.size(); i++) {
		Eigen::Vector3d vert = m.V.row(boundary_vertices[i]);
		double d = (vert - start_v).norm();
		if (d < min) {
			min = d;
			index = i;
		}
	}

	vector<int> reordered;
	for (int i = 0; i < boundary_vertices.size(); i++) {
		reordered.push_back(boundary_vertices[(i + index) % boundary_vertices.size()]);
	}
	return reordered;
}

/*void LaplacianRemesh::reverse_path(vector<int> path_vertices) {
	reverse(path_vertices.begin(), path_vertices.end());
	int v = path_vertices[path_vertices.size() - 1];
	path_vertices.insert(path_vertices.begin() + 1, path_vertices.begin(), path_vertices.end() - 1);
	path_vertices[0] = v;
}*/

/** For cut we determine if the boundary is CCW by comparing the direction of the normal of the remaining boundary points against the vector from the center of the removed boundary points to the center of the remaining boundary points.
	For extrusion we simply look at the sign of the inner area of the projected points (since we're looking from the "correct" viewpoint)
**/
bool LaplacianRemesh::is_counter_clockwise_boundaries(Eigen::MatrixXd boundary_points, Eigen::Matrix4f modelview, Eigen::Matrix4f proj, Eigen::Vector4f viewport, Eigen::RowVector3d mean_viewpoint, bool cut) {
	if (cut) {
		Eigen::RowVector3d center = boundary_points.colwise().mean();
		Eigen::Vector3d normal(0, 0, 0);
		Eigen::Vector3d vec0, vec1;
		for (int i = 0; i < boundary_points.rows(); i++) {
			vec0 = boundary_points.row(i) - center;
			vec1 = boundary_points.row((i + 1) % boundary_points.rows()) - center;
			normal += vec1.cross(vec0);
		}
		normal.normalize();

		if (normal.dot(center - mean_viewpoint) < 0) {
			return false;
		}
		else {
			return true;
		}
	}
	else {
		double total_area = 0.0;
		Eigen::RowVector3d pt, vert;
		Eigen::RowVector2d prev, next;
		vert = boundary_points.row(boundary_points.rows() - 1);
		igl::project(vert, modelview, proj, viewport, pt);
		prev = pt.leftCols(2);

		for (int i = 0; i < boundary_points.rows(); i++) {
			vert = boundary_points.row(i);
			igl::project(vert, modelview, proj, viewport, pt);
			next = pt.leftCols(2);
			total_area += (prev[1] + next[1]) * (next[0] - prev[0]);
			prev = next;
		}

		if (total_area > 0) {
			return false;
		}

		return true;
	}
}

double LaplacianRemesh::compute_average_distance_between_onPolygon_vertices(std::vector<PathElement> path) {
	Eigen::MatrixX3d onPoly_vertices = Eigen::MatrixX3d::Zero(0, 3);
	for (int i = 0; i < path.size(); i++) {
		if (path[i].get_type() == PathElement::ElementType::FACE) {
			onPoly_vertices.conservativeResize(onPoly_vertices.rows() + 1, Eigen::NoChange);
			onPoly_vertices.row(onPoly_vertices.rows() - 1) = path[i].get_vertex().transpose();
		}
	}

	double total = 0.0;
	for (int i = 0; i < onPoly_vertices.rows() - 1; i++) {
		total += (onPoly_vertices.row(i) - onPoly_vertices.row(i + 1)).norm();
	}
	return total / (onPoly_vertices.rows() - 1);
}

double LaplacianRemesh::compute_average_length_of_crossing_edges(std::vector<PathElement> path, Eigen::MatrixXd startV, Eigen::MatrixXi startEV) {
	double total = 0.0;
	int count = 0;
	for (int i = 0; i < path.size(); i++) {
		if (path[i].get_type() == PathElement::ElementType::EDGE) {
			total += (startV.row(startEV(path[i].get_ID(), 0)) - startV.row(startEV(path[i].get_ID(), 1))).norm();
			count++;
		}
	}
	return total / count;
}


void LaplacianRemesh::update_open_path_points_and_bindings(Stroke& stroke, SurfacePath& surface_path) {
	vector<PathElement> path = surface_path.get_path();
	int nr_to_remove = 0;
	if (path[0].get_vertex() == path[path.size() - 1].get_vertex()) {
		nr_to_remove = 1;
	}

	//TODO: check if we actually want a loop below
	Eigen::MatrixX3d new_3DPoints(path.size() + 1 - nr_to_remove, 3); //Increase size by 1 (if it's not a loop) because we want it to become a loop again
	vector<int> new_closest_vertex_indices(path.size() + 1 - nr_to_remove);
	for (int i = 0; i < path.size() - nr_to_remove; i++) {
		new_3DPoints.row(i) = path[i].get_vertex().transpose();
		new_closest_vertex_indices[i] = path[i].get_v_idx();
	}
	new_3DPoints.row(new_3DPoints.rows() - 1) = new_3DPoints.row(0);
	new_closest_vertex_indices[new_closest_vertex_indices.size() - 1] = new_closest_vertex_indices[0];

	stroke.set3DPoints(new_3DPoints);
	stroke.set_closest_vert_bindings(new_closest_vertex_indices);
}
