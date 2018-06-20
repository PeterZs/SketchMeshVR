#include "Patch.h"
#include <igl/edge_topology.h>

static Eigen::MatrixXi EV, FE, EF;
static int next_mesh_ID;

/*Patch& Patch::operator=(Patch other) {
	mesh = other.mesh;
	parent_vertices = other.parent_vertices;
	parent_faces = other.parent_faces;
	new_vertex_boundary_markers = other.new_vertex_boundary_markers;
	new_part_of_original_stroke = other.new_part_of_original_stroke;
	//new_mapped_indices = other.new_mapped_indices;
	return *this;
}*/

std::vector<Patch*> Patch::init_patches(Mesh& h) {
	next_mesh_ID = h.ID + 1;
	igl::edge_topology(h.V, h.F, EV, FE, EF);
	
	std::vector<Patch*> patches;
	h.face_patch_map.resize(h.F.rows(), nullptr);
	
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
		int face2 = (EF(edge, 0) == face) ? EF(edge, 1) : EF(edge, 0); //get the adjacent polygon
		if (face_patch_map[face2] == nullptr && !sharp_edge[edge]) {
			propagate_patch(patch, face2, faces, face_patch_map, sharp_edge);
		}
	}
}

void Patch::create_mesh_structure(Mesh& m, Eigen::VectorXi& faces) {
	
	patch_vertices.resize(0, 3);
	patch_vertex_is_init = Eigen::VectorXi::Zero(m.V.rows());
	patch_faces.resize(0,3);
	patch_edge_is_init = Eigen::VectorXi::Zero(EF.rows());
	for (int i = 0; i < faces.rows(); i++) {
		for (int j = 0; j < 3; j++) {
			int edge = FE(i, j);
			get_patch_edge(edge, patch_edge_is_init, patch_vertex_is_init, i, patch_vertices, m.sharp_edge, new_sharp_edge, m.V, m.vertex_boundary_markers, m.part_of_original_stroke, m.new_mapped_indices);
		}
		patch_faces.conservativeResize(patch_faces.rows() + 1, Eigen::NoChange);
		patch_faces.bottomRows(1) << m.F.row(faces[i]);
		parent_faces.conservativeResize(parent_faces.rows() + 1, Eigen::NoChange);
		parent_faces[i] = faces[i];	
	}
	mesh.V = patch_vertices;
	mesh.F = patch_faces;
	mesh.vertex_boundary_markers = new_vertex_boundary_markers;
	mesh.part_of_original_stroke = new_part_of_original_stroke;
	mesh.new_mapped_indices = new_mapped_indices;
	mesh.sharp_edge = new_sharp_edge;
	mesh.ID = next_mesh_ID;

	next_mesh_ID++;
}

void Patch::get_patch_edge(int edge, Eigen::VectorXi& patch_edge_is_init, Eigen::VectorXi& patch_vertex_is_init, int face, Eigen::MatrixXd& patch_vertices, Eigen::VectorXi& sharp_edge, Eigen::VectorXi& new_sharp_edge, Eigen::MatrixXd& V_orig, Eigen::VectorXi& boundary_markers_orig, Eigen::VectorXi& part_of_original_orig, Eigen::VectorXi& new_mapped_indices_orig) {
	if (sharp_edge(edge)) {
		//Do nothing
	}
	else if (patch_edge_is_init[edge]) {
		return;
	}

	int start = EV(edge, 0);
	int end = EV(edge, 1);
	get_patch_vertex(start, face, patch_vertices, patch_vertex_is_init, V_orig, boundary_markers_orig, part_of_original_orig, new_mapped_indices_orig);
	get_patch_vertex(end, face, patch_vertices, patch_vertex_is_init, V_orig, boundary_markers_orig, part_of_original_orig, new_mapped_indices_orig);
	new_sharp_edge.conservativeResize(new_sharp_edge.rows() + 1);
	new_sharp_edge.tail(1) << sharp_edge[edge];
	patch_edge_is_init[edge] = 1;

}


void Patch::get_patch_vertex(int v_idx, int face, Eigen::MatrixXd& patch_vertices, Eigen::VectorXi& patch_vertex_is_init, Eigen::MatrixXd& V_orig, Eigen::VectorXi& boundary_markers_orig, Eigen::VectorXi& part_of_original_orig, Eigen::VectorXi& new_mapped_indices_orig) {
	if (patch_vertex_is_init[v_idx]) {
		return;
	}

	patch_vertices.conservativeResize(patch_vertices.rows() + 1, Eigen::NoChange);
	patch_vertices.bottomRows(1) = V_orig.row(v_idx);
	parent_vertices.conservativeResize(parent_vertices.rows() + 1);
	parent_vertices.tail(1) << v_idx;
	new_vertex_boundary_markers.conservativeResize(new_vertex_boundary_markers.rows() + 1);
	new_vertex_boundary_markers.tail(1) << boundary_markers_orig.row(v_idx);
	new_part_of_original_stroke.conservativeResize(new_part_of_original_stroke.rows() + 1);
	new_part_of_original_stroke.tail(1) << part_of_original_orig.row(v_idx);
	new_mapped_indices.conservativeResize(new_mapped_indices.rows() + 1);
	new_mapped_indices.tail(1) << -1;// new_mapped_indices_orig.row(v_idx); //TODO: CHECK IF WE NEED THE NEW_MAPPED_INDICES OR IF THESE ARE ALWAYS UNINITIALIZED WHEN THIS IS CALLED
	patch_vertex_is_init[v_idx] = 1;
}

void Patch::update_parent_vertex_positions() {
	for (int i = 0; i < mesh.V.rows(); i++) {
		V[parent_vertices[i]] = mesh.V.row(i);
	}
}

/*
public Polyhedron mesh;
// split at internal sharp edges
public void create_mesh_structure(Polyhedron h, ArrayList faces) {
	mesh = new Polyhedron();	

	PatchEdge[] patch_edges = new PatchEdge[h.edges.size()];
	PatchVertex[] patch_vertices = new PatchVertex[h.faces.size() * 3];

	for (int i = 0; i<faces.size(); i++) {
		Face face = (Face)faces.get(i);

		PatchEdge[] edges = new PatchEdge[3];
		for (int j = 0; j<3; j++) {
			Edge edge = face.edges[j];
			edges[j] = get_patch_edge(edge, patch_edges, face, patch_vertices);
		}
		PatchFace patchface = new PatchFace(edges[0], edges[1], edges[2]);
		patchface.parent_face = face;
		mesh.append(patchface);
	}
}

public PatchEdge get_patch_edge(Edge edge, PatchEdge[] patch_edges, Face face, PatchVertex[] patch_vertices) {
	if (edge.sharp || edge.get_another_polygon(face) == null) {
	}
	else if (patch_edges[edge.index] != null)
		return patch_edges[edge.index];

	PatchVertex start = get_patch_vertex(edge.start, face, patch_vertices);
	PatchVertex end = get_patch_vertex(edge.end, face, patch_vertices);

	PatchEdge patch_edge = new PatchEdge(start, end);
	patch_edge.parent_edge = edge;
	patch_edge.seam = edge.seam;
	patch_edges[edge.index] = patch_edge;

	mesh.append(patch_edge);

	return patch_edge;

}
public PatchVertex get_patch_vertex(Vertex v, Face face, PatchVertex[] patch_vertices) {
	if (patch_vertices[face.index * 3 + face.get_vertex_index(v)] != null)
		return patch_vertices[face.index * 3 + face.get_vertex_index(v)];

	PatchVertex pv = new PatchVertex(v.x, v.y, v.z);
	pv.normal = v.normal;
	pv.parent_vertex = v;
	pv.fixed = v.fixed;
	pv.prescribed_normal = v.prescribed_normal;

	mesh.append(pv);

	// fill
	ArrayList pairs = get_vertex_face_pairs(v, face);
	for (int i = 0; i<pairs.size(); i++) {
		Object[] pair = (Object[]) pairs.get(i);
		Vertex u = (Vertex)pair[0];
		Face f = (Face)pair[1];
		patch_vertices[f.index * 3 + f.get_vertex_index(u)] = pv;
	}

	return pv;
}

// fan propagate
public ArrayList get_vertex_face_pairs(Vertex v, Face face) {

	ArrayList pairs = new ArrayList();

	if (!v.on_sharp()) {
		ArrayList faces = v.polygons();
		for (int i = 0; i<faces.size(); i++) {
			Object[] pair = { v, faces.get(i) };
			pairs.add(pair);
		}
		return pairs;
	}


	get_vertex_face_pairs_sub(v, face, null, pairs);
	//get_sharp_normal_sub(polygon.edges[i], polygon, polygons, polygon);   

	return pairs;
}
public void get_vertex_face_pairs_sub(Vertex v, Face face, Edge prev_edge, ArrayList pairs) {
	Object[] pair = { v, face };
	pairs.add(pair);

	for (int i = 0; i<3; i++) {
		Edge next_edge = face.edges[i];
		if (next_edge.contain(v) && next_edge != prev_edge && !next_edge.sharp) {
			Face next_face = next_edge.get_another_polygon(face);
			if (next_face != null)
				get_vertex_face_pairs_sub(v, next_face, next_edge, pairs);
		}
	}
}
*/