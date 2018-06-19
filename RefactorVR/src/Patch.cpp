#include "Patch.h"


static std::vector<Patch*> init_patches(Mesh& h) {
	std::vector<Patch*> patches(h.F.rows());
	for (int i = 0; i < h.F.rows(); i++) {
		Patch* face_patch = h.face_patch_map[i];
		if (face_patch == nullptr) {

			patches.push_back(face_patch);
		}
	}
	return patches;
}

public static ArrayList init_patches(Polyhedron h) {
	h.set_indices();

	ArrayList patches = new ArrayList();
	for (int i = 0; i<h.faces.size(); i++) {
		h.faces(i).patch = null;
	}
	for (int i = 0; i<h.faces.size(); i++) {
		Face face = h.faces(i);
		if (face.patch == null) {
			Patch patch = new Patch();
			ArrayList faces = new ArrayList();
			propagate_patch(patch, face, faces);
			patch.create_mesh_structure(h, faces);
			patches.add(patch);
		}
	}

	return patches;
}



static void propagate_patch(Patch patch, Face face, ArrayList faces) {
	face.patch = patch;
	faces.add(face);
	for (int i = 0; i<3; i++) {
		Edge edge = face.edges[i];
		Face face2 = edge.get_the_other_polygon(face);
		if (face2 != null && face2.patch == null && !edge.sharp) {
			propagate_patch(patch, face2, faces);
		}
	}
}




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
