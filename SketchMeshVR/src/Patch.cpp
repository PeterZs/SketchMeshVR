#include "Patch.h"
#include <igl/edge_topology.h>
#include <igl/internal_angles.h>
#include <igl/upsample.h>
#include <igl/slice.h>
#include <igl/adjacency_list.h>
#include <igl/cat.h>
#include <unordered_set>

static Eigen::MatrixXi EV, FE, EF, patchEV, patchFE, patchEF;
static int next_mesh_ID;

std::vector<Patch*> Patch::init_patches(Mesh& h) {
	next_mesh_ID = h.ID + 1;
	igl::edge_topology(h.V, h.F, EV, FE, EF);

	std::vector<Patch*> patches;
	h.face_patch_map.clear();
	h.face_patch_map.resize(h.F.rows(), nullptr);

	std::cout << "Rows:" << h.sharp_edge.rows() << std::endl;;

	for (int i = 0; i < h.F.rows(); i++) {
		Patch* face_patch = h.face_patch_map[i];
		if (face_patch == nullptr) {
			Patch* new_patch = new Patch();
			Eigen::VectorXi faces;
			propagate_patch(new_patch, i, faces, h.face_patch_map, h.sharp_edge);
			(*new_patch).create_mesh_structure(h, faces);
			patches.push_back(new_patch);
		}
	}
	return patches;
}

void Patch::propagate_patch(Patch* patch, int face, Eigen::VectorXi& faces, std::vector<Patch*> &face_patch_map, Eigen::VectorXi& sharp_edge) {
	face_patch_map[face] = patch;
	faces.conservativeResize(faces.rows() + 1);
	faces.tail(1) << face;
	for (int i = 0; i < 3; i++) {
		int edge = FE(face, i);
		int face2 = (EF(edge, 0) == face) ? EF(edge, 1) : EF(edge, 0); //get the adjacent polygon.
		if (face2 != -1 && (face_patch_map[face2] == nullptr) && (!sharp_edge[edge])) {
			propagate_patch(patch, face2, faces, face_patch_map, sharp_edge);
		}
	}
}

void Patch::create_mesh_structure(Mesh& m, Eigen::VectorXi& faces) {
	patch_vertices.resize(0, 3);
	patch_vertex_is_init = Eigen::VectorXi::Zero(m.V.rows());
	patch_faces.resize(0, 3);
	patch_edge_is_init = Eigen::VectorXi::Zero(EF.rows());
	mesh_to_patch_indices = Eigen::VectorXi::Constant(m.V.rows(), -1);
	Eigen::MatrixXi tmp_sharp_edge(0, 2), tmp_edge_boundary_markers(0, 3);

	for (int i = 0; i < faces.rows(); i++) {
		for (int j = 0; j < 3; j++) {
			int edge = FE(faces[i], j);
			get_patch_edge(edge, patch_edge_is_init, patch_vertex_is_init, faces[i], patch_vertices, m.sharp_edge, tmp_sharp_edge, m.edge_boundary_markers, tmp_edge_boundary_markers, m.V, m.vertex_is_fixed, m.new_mapped_indices, mesh_to_patch_indices);
		}
	}

	for (int i = 0; i < faces.rows(); i++) {
		patch_faces.conservativeResize(patch_faces.rows() + 1, Eigen::NoChange);
		patch_faces.bottomRows(1) << mesh_to_patch_indices[m.F(faces[i], 0)], mesh_to_patch_indices[m.F(faces[i], 1)], mesh_to_patch_indices[m.F(faces[i], 2)];
	}

	igl::edge_topology(patch_vertices, patch_faces, patchEV, patchFE, patchEF);

	Eigen::VectorXi col1Equals, col2Equals;
	int start, end, val, equal_pos;
	new_sharp_edge.resize(patchEV.rows());
	new_sharp_edge.setZero();
	for (int i = 0; i < tmp_sharp_edge.rows(); i++) { //Go over all sharp edges in the patch and set them
		start = mesh_to_patch_indices[tmp_sharp_edge(i, 0)];
		end = mesh_to_patch_indices[tmp_sharp_edge(i, 1)];
		col1Equals = patchEV.col(0).cwiseEqual(std::min(start, end)).cast<int>();
		col2Equals = patchEV.col(1).cwiseEqual(std::max(start, end)).cast<int>();
		val = (col1Equals + col2Equals).maxCoeff(&equal_pos);
		if (val == 2) {
			new_sharp_edge[equal_pos] = 1;
		}
	}

	new_edge_boundary_markers.resize(patchEV.rows());
	for (int i = 0; i < tmp_edge_boundary_markers.rows(); i++) {
		start = mesh_to_patch_indices[tmp_edge_boundary_markers(i, 0)];
		end = mesh_to_patch_indices[tmp_edge_boundary_markers(i, 1)];
		col1Equals = patchEV.col(0).cwiseEqual(std::min(start, end)).cast<int>();
		col2Equals = patchEV.col(1).cwiseEqual(std::max(start, end)).cast<int>();
		val = (col1Equals + col2Equals).maxCoeff(&equal_pos);
		if (val == 2) {
			new_edge_boundary_markers[equal_pos] = tmp_edge_boundary_markers(i, 2);
		}
	}

	mesh.V = patch_vertices;
	mesh.F = patch_faces;
	mesh.edge_boundary_markers = new_edge_boundary_markers;
	mesh.vertex_is_fixed = new_vertex_is_fixed;
	mesh.new_mapped_indices = new_mapped_indices;
	mesh.sharp_edge = new_sharp_edge;
	mesh.ID = next_mesh_ID;
	mesh.mesh_to_patch_indices = mesh_to_patch_indices; //Maps a vertex index from the base mesh to a vertex index in the mesh of a underlying patch
	next_mesh_ID++;
}

void Patch::get_patch_edge(int edge, Eigen::VectorXi& patch_edge_is_init, Eigen::VectorXi& patch_vertex_is_init, int face, Eigen::MatrixXd& patch_vertices, Eigen::VectorXi& sharp_edge, Eigen::MatrixXi& tmp_sharp_edge, Eigen::VectorXi &edge_boundary_markers, Eigen::MatrixXi &tmp_edge_boundary_markers, Eigen::MatrixXd& V_orig, Eigen::VectorXi & vertex_is_fixed_orig, Eigen::VectorXi& new_mapped_indices_orig, Eigen::VectorXi& mesh_to_patch_indices) {
	if (sharp_edge(edge)) {

	}
	if (patch_edge_is_init[edge]) {
		return;
	}

	int start = EV(edge, 0);
	int end = EV(edge, 1);
	get_patch_vertex(start, face, patch_vertices, patch_vertex_is_init, V_orig, vertex_is_fixed_orig, new_mapped_indices_orig, mesh_to_patch_indices);
	get_patch_vertex(end, face, patch_vertices, patch_vertex_is_init, V_orig, vertex_is_fixed_orig, new_mapped_indices_orig, mesh_to_patch_indices);
	if (sharp_edge[edge]) {
		tmp_sharp_edge.conservativeResize(tmp_sharp_edge.rows() + 1, Eigen::NoChange);
		tmp_sharp_edge.bottomRows(1) << start, end;
	}
	tmp_edge_boundary_markers.conservativeResize(tmp_edge_boundary_markers.rows() + 1, Eigen::NoChange);
	tmp_edge_boundary_markers.bottomRows(1) << start, end, edge_boundary_markers[edge];
	patch_edge_is_init[edge] = 1;
}

void Patch::get_patch_vertex(int v_idx, int face, Eigen::MatrixXd& patch_vertices, Eigen::VectorXi& patch_vertex_is_init, Eigen::MatrixXd& V_orig, Eigen::VectorXi& vertex_is_fixed_orig, Eigen::VectorXi& new_mapped_indices_orig, Eigen::VectorXi& mesh_to_patch_indices) {
	if (patch_vertex_is_init[v_idx]) {
		return;
	}

	patch_vertices.conservativeResize(patch_vertices.rows() + 1, Eigen::NoChange);
	patch_vertices.bottomRows(1) = V_orig.row(v_idx);
	parent_vertices.conservativeResize(parent_vertices.rows() + 1);
	parent_vertices.tail(1) << v_idx;
	mesh_to_patch_indices[v_idx] = parent_vertices.rows() - 1;
	new_vertex_is_fixed.conservativeResize(new_vertex_is_fixed.rows() + 1);
	new_vertex_is_fixed.tail(1) << vertex_is_fixed_orig.row(v_idx);
	new_mapped_indices.conservativeResize(new_mapped_indices.rows() + 1);
	new_mapped_indices.tail(1) << -1;
	patch_vertex_is_init[v_idx] = 1;
}

void Patch::update_parent_vertex_positions(Eigen::MatrixXd& base_V) {
	for (int i = 0; i < mesh.V.rows(); i++) {
		base_V.row(parent_vertices[i]) = mesh.V.row(i);
	}
}

void Patch::update_patch_vertex_positions(Eigen::MatrixXd& base_V) {
	for (int i = 0; i < base_V.rows(); i++) {
		if (mesh_to_patch_indices[i] != -1) {
			mesh.V.row(mesh_to_patch_indices[i]) = base_V.row(i);
		}
	}
}

double Patch::compute_mean_minimum_triangle_angle() {
	Eigen::MatrixXd triangle_angles;
	igl::internal_angles(mesh.V, mesh.F, triangle_angles);
	Eigen::VectorXd min_angles_per_triangle = triangle_angles.rowwise().minCoeff();
	return min_angles_per_triangle.mean();
}

void Patch::upsample_patch(Eigen::MatrixXd& base_V, Eigen::MatrixXi& base_F, std::vector<Patch*> base_face_patch_map, Eigen::VectorXi& base_edge_boundary_markers, Eigen::VectorXi& base_sharp_edge, Eigen::VectorXi& base_vertex_is_fixed, Eigen::MatrixXi& replacing_vertex_bindings) {
	replacing_vertex_bindings.resize(0, 4);
	Eigen::MatrixXi old_baseF = base_F;
	int start_size_baseV = base_V.rows();
	//Store EV of base_mesh before upsampling
	Eigen::MatrixXi startEV, startFE, startEF;
	igl::edge_topology(base_V, base_F, startEV, startFE, startEF);

	Eigen::MatrixXd new_patchV;
	Eigen::MatrixXi new_patchF;
	igl::upsample(mesh.V, mesh.F, new_patchV, new_patchF);
	//Assuming that each face is only part of 1 Patch, we can take the base_mesh's face_patch_map, and remove slice (keeping only the faces that do not map to the current patch)
	std::vector<int> clean_faces;
	for (int i = 0; i < base_face_patch_map.size(); i++) {
		if (base_face_patch_map[i] == this) {
			//Face should be removed from base_F		
		}
		else {
			clean_faces.push_back(i);
		}
	}

	//Store all clean faces from the base mesh (that are not part of this Patch) in tmp_F
	Eigen::MatrixXi tmp_F;
	Eigen::VectorXi row_idx = Eigen::VectorXi::Map(clean_faces.data(), clean_faces.size());
	Eigen::VectorXi col_idx(3);
	col_idx.col(0) << 0, 1, 2;
	igl::slice(base_F, row_idx, col_idx, tmp_F);

	//Append newly created vertices to the vertex list of the base (don't care about the patch, because this will be reinitialized afterwards anyway)
	int nr_added_points = new_patchV.rows() - mesh.V.rows();
	Eigen::MatrixXd newly_addedV = new_patchV.bottomRows(nr_added_points);
	base_V.conservativeResize(base_V.rows() + nr_added_points, Eigen::NoChange);
	base_V.bottomRows(nr_added_points) = newly_addedV;

	int start_size = mesh.V.rows();
	parent_vertices.conservativeResize(parent_vertices.rows() + nr_added_points);
	for (int i = 0; i < nr_added_points; i++) {
		parent_vertices[start_size + i] = start_size_baseV + i;
	}

	Eigen::MatrixXi new_patchF_reindexed(new_patchF.rows(), 3);
	for (int i = 0; i < new_patchF.rows(); i++) {
		new_patchF_reindexed.row(i) << parent_vertices[new_patchF(i, 0)], parent_vertices[new_patchF(i, 1)], parent_vertices[new_patchF(i, 2)];
	}

	//Append newly created faces to the cleaned-up face list of the base
	base_F.resize(tmp_F.rows() + new_patchF_reindexed.rows(), 3);
	base_F.topRows(tmp_F.rows()) = tmp_F;
	base_F.bottomRows(new_patchF_reindexed.rows()) = new_patchF_reindexed;


	//Construct EV for the new base_V and base_F
	//Loop over all these edges. If there was a connection between the 2 edge vertices before upsampling, copy the edge_boundary_marker for the specific edge. 
	//Else it is a subdivided edge. If exactly one of the edge vertices already existed in the start_V, then we have an edge that is half of a broken up edge (if both vertices are new, it is an internal edge)
	//For the subdivided edge, find the other subdivided edge that attaches to the middle vertex. Check in the original vertex_boundary_markers if the edge between these 2 endpoints was marked, and copy that marker onto both edges

	Eigen::MatrixXi newEV, newFE, newEF;
	igl::edge_topology(base_V, base_F, newEV, newFE, newEF);

	//Remove faces that are on adjacent patches, and are attached to edges that will be split in the current patch
	//Eigen::MatrixXi has_missing_face = newEF.cwiseEqual(-1).cast<int>();
	//Eigen::VectorXi has_missing = (has_missing_face.rowwise().sum().array() > 0).cast<int>();

	/*int present_face;
	for (int i = 0; i < newEF.rows(); i++) {
		if (has_missing[i]) {
			present_face = (newEF(i, 0) == -1 ? newEF(i, 1) : newEF(i, 0));
			if (base_F(present_face, 0) < start_size_baseV && base_F(present_face, 1) < start_size_baseV && base_F(present_face, 2) < start_size_baseV) {
				dirty_faces[present_face] = 1;
			}
		}
	}*/


	//Update topology again with the previous patch-adjacent faces (with split edges) removed
	//igl::edge_topology(base_V, base_F, newEV, newFE, newEF);



	Eigen::VectorXi new_base_edge_boundary_markers(newEV.rows()), new_base_sharp_edge(newEV.rows());
	new_base_edge_boundary_markers.setConstant(-1);
	new_base_sharp_edge.setConstant(-1);
	Eigen::VectorXi new_base_vertex_is_fixed(base_vertex_is_fixed.rows() + nr_added_points);
	new_base_vertex_is_fixed.topRows(base_vertex_is_fixed.rows()) = base_vertex_is_fixed;

	//std::vector<bool> dirty_face(base_F.rows()); //bools are defaulted to false

	Eigen::MatrixXi faces_to_add(0, 3);
	Eigen::VectorXi dirty_faces(base_F.rows());
	dirty_faces.setZero();
	int dirty_count = 0;
	for (int i = 0; i < newEV.rows(); i++) {
		if (new_base_sharp_edge[i] != -1) { //Edge has already been covered
			continue;
		}

		int start = newEV(i, 0);
		int end = newEV(i, 1);
		if (start >= start_size_baseV && end >= start_size_baseV) { //Both vertices are newly added, so they form an interior edge.
			new_base_edge_boundary_markers[i] = 0;
			new_base_sharp_edge[i] = 0;
		}
		else if (start < start_size_baseV && end < start_size_baseV) { //Both vertices were existing already
			std::cout << "This shouldn't happen (yet) " << start << "  " << end << std::endl;
			int old_edge = find_edge(start, end, startEV);
			new_base_edge_boundary_markers[i] = base_edge_boundary_markers[old_edge];
			new_base_sharp_edge[i] = base_sharp_edge[old_edge];
		}
		else { //One vertex is old (smallest index), one is new (larger index)
			int equal_pos;
			//Find all edges in newEV that connect to the vertex with the larger index (new vertex). Only one of these (except for the current edge i) will be with a vertex that has an index < start_size_baseV
			Eigen::VectorXi right_end_point = (newEV.col(1).cwiseEqual(std::max(start, end))).cast<int>();
			Eigen::VectorXi connects_old_point = (newEV.col(0).array() < start_size_baseV).cast<int>();
			Eigen::VectorXi not_current_edge = (newEV.col(0).array() != std::min(start, end)).cast<int>();
			int maxval = (right_end_point + connects_old_point + not_current_edge).maxCoeff(&equal_pos); //Find the row that contains both vertices of this edge
			if (maxval == 3) {
				int old_edge_start = std::min(start, end);
				int old_edge_end = std::min(newEV(equal_pos, 0), newEV(equal_pos, 1));
				int original_edge = find_edge(old_edge_start, old_edge_end, startEV);
				int other_edge = find_edge(old_edge_end, std::max(start, end), newEV);
				int marker_val = base_edge_boundary_markers[original_edge];
				int sharp_val = base_sharp_edge[original_edge];
				new_base_edge_boundary_markers[i] = marker_val;
				new_base_sharp_edge[i] = sharp_val;
				new_base_edge_boundary_markers[other_edge] = marker_val;
				new_base_sharp_edge[other_edge] = sharp_val;


				int tmp1 = std::min(old_edge_start, old_edge_end);
				int tmp2 = std::max(old_edge_start, old_edge_end);

				replacing_vertex_bindings.conservativeResize(replacing_vertex_bindings.rows() + 1, Eigen::NoChange);
				replacing_vertex_bindings.bottomRows(1) << marker_val, tmp1, tmp2, std::max(start, end);

				//Only set the newly added midpoint vertices to fixed if both endpoints of the original edge are fixed and the edge has a boundary marker (to prevent fixed edges from forming when 2 endpoints from the original edge each are fixed but are not connected by a boundary marked edge)
				if (base_vertex_is_fixed[old_edge_start] && base_vertex_is_fixed[old_edge_end] && base_edge_boundary_markers[original_edge]) {
					new_base_vertex_is_fixed[std::max(start, end)] = 1;
				}
				else {
					new_base_vertex_is_fixed[std::max(start, end)] = 0;
				}


				//Find face that is adjacent to the patch boundary edge that we just split, remove that face and replace it with 2 half faces (outside the patch's domain)
				if (sharp_val) {
					int edge = find_edge(old_edge_start, old_edge_end, startEV);
					int cur_edge = find_edge(old_edge_start, old_edge_end, newEV); //Find the edge between the original triangle "base" in the current EV
					int face1 = startEF(edge, 0);
					int face2 = startEF(edge, 1);

					int mid_vertex = std::max(start, end);
					int vert_to_connect;
					if (base_face_patch_map[face1] == this) { //Patch face that is removed, take the other one
						for (int j = 0; j < 3; j++) {
							if ((old_baseF(face2, j) != old_edge_start) && (old_baseF(face2, j) != old_edge_end)) {
								vert_to_connect = old_baseF(face2, j);
								faces_to_add.conservativeResize(faces_to_add.rows() + 1, Eigen::NoChange);
								faces_to_add.bottomRows(1) << old_edge_start, vert_to_connect, mid_vertex;
								faces_to_add.conservativeResize(faces_to_add.rows() + 1, Eigen::NoChange);
								faces_to_add.bottomRows(1) << mid_vertex, vert_to_connect, old_edge_end;
								int face2_cur = (newEF(cur_edge, 0) == -1) ? newEF(cur_edge, 1) : newEF(cur_edge, 0);
								std::cout << newEF.row(cur_edge) << "   " << base_F.row(face2_cur) << std::endl;
								dirty_faces[face2_cur] = 1;
							//	dirty_count++;
							}
						}
					}
					else if (base_face_patch_map[face2] == this) {
						for (int j = 0; j < 3; j++) {
							if ((old_baseF(face1, j) != old_edge_start) && (old_baseF(face1, j) != old_edge_end)) {
								vert_to_connect = old_baseF(face1, j);
								faces_to_add.conservativeResize(faces_to_add.rows() + 1, Eigen::NoChange);
								faces_to_add.bottomRows(1) << old_edge_start, mid_vertex, vert_to_connect;
								faces_to_add.conservativeResize(faces_to_add.rows() + 1, Eigen::NoChange);
								faces_to_add.bottomRows(1) << mid_vertex, old_edge_end, vert_to_connect;
								int face1_cur = (newEF(cur_edge, 0) == -1) ? newEF(cur_edge, 1) : newEF(cur_edge, 0);
								dirty_faces[face1_cur] = 1;
								//dirty_count++;
							}
						}
					}
					else {
						std::cerr << "Something went wrong? " << std::endl;
					}
					std::cout << "Test: " << old_edge_start << "  " << mid_vertex << "   " << vert_to_connect << "   " << old_edge_end << std::endl;
				}
				
			}else {
				std::cout << "What's up: " << start << "  " << end << std::endl;
			}
		}
	}

	//std::cout << "Dirty count: " << dirty_count << std::endl;

	clean_faces.clear();
	for (int i = 0; i < base_F.rows(); i++) {
		if (!dirty_faces[i]) {
			clean_faces.push_back(i);
		}
	}

	row_idx = Eigen::VectorXi::Map(clean_faces.data(), clean_faces.size());
	igl::slice(base_F, row_idx, col_idx, tmp_F);
	base_F = tmp_F;

	base_F = igl::cat(1, base_F, faces_to_add);

	base_sharp_edge = new_base_sharp_edge;
	base_edge_boundary_markers = new_base_edge_boundary_markers;
	base_vertex_is_fixed = new_base_vertex_is_fixed;
}

int Patch::find_edge(int start, int end, Eigen::MatrixXi& EV) {
	Eigen::VectorXi col1Equals, col2Equals;
	int equal_pos;
	col1Equals = EV.col(0).cwiseEqual(std::min(start, end)).cast<int>();
	col2Equals = EV.col(1).cwiseEqual(std::max(start, end)).cast<int>();
	int val = (col1Equals + col2Equals).maxCoeff(&equal_pos); //Find the row that contains both vertices of this edge

	if (val != 2) {
		return -1; //Edge between start and end doesn't exist in EV
	}
	return equal_pos;
}