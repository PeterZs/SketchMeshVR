#include <igl/unproject_onto_mesh.h>
#include <igl/edge_topology.h>
#include "SurfacePath.h"
#include "Plane.h"

Eigen::MatrixXi EV, FE, EF;
Eigen::MatrixXd N_faces;
Stroke& origin_stroke;

using namespace std;

SurfacePath::SurfacePath() {

}

//Adds stroke elements at the intersection points of the original drawn stroke with mesh edges
void SurfacePath::create_from_stroke(Stroke & stroke) {
	origin_stroke = stroke;
	Eigen::Matrix4f modelview = stroke.viewer.core.view * stroke.viewer.core.model;
	int faceID = -1;
	Eigen::Vector3f bc;

	int prev_p = 0;
	int start_p = prev_p;
	int next_p;
	igl::unproject_onto_mesh(stroke.get_stroke2DPoints().row(prev_p), modelview, stroke.viewer.core.proj, stroke.viewer.core.viewport, stroke.get_V(), stroke.get_F(), faceID, bc);

	int start_face = faceID;
	int n = 1;
	Eigen::RowVector3d pt(0, 0, 0);

	igl::edge_topology(stroke.get_V(), stroke.get_F(), EV, FE, EF);
	igl::per_face_normals(stroke.get_V(), stroke.get_F(), N_faces);

	while(true) {
		next_p = n;
		pt = stroke.get_V().row(stroke.get_F()(faceID, 0))*bc(0) + stroke.get_V().row(stroke.get_F()(faceID, 1))*bc(1) + stroke.get_V().row(stroke.get_F()(faceID, 2))*bc(2);
		PathElement newElement(faceID, PathElement::FACE, pt);
		path.push_back(newElement);

		bool forward;
		faceID = extend_path(prev_p, next_p, faceID, forward, false, modelview);

		if(!forward) {
			next_p = prev_p;
		}

		if(next_p == start_p && faceID == start_face) {
			break;
		}

		if(front_facing(faceID)) {
			n = next_p + 1;
		} else {
			n = next_p - 1;
		}
		prev_p = next_p;
	}
}

int SurfacePath::extend_path(int prev_p, int next_p, int faceID, bool forward, bool front_facing_only, Eigen::Matrix4f modelview) {
	Eigen::Vector3d source, dir;
	Eigen::MatrixX2d stroke2DPoints = origin_stroke.get_stroke2DPoints();
	Eigen::MatrixX3d stroke3DPoints = origin_stroke.get3DPoints();
	Eigen::Vector2d tmp = stroke2DPoints.row(prev_p);
	igl::unproject_ray(tmp, modelview, origin_stroke.viewer.core.proj, origin_stroke.viewer.core.viewport, source, dir);
	Plane cutPlane(source, stroke3DPoints.row(prev_p), stroke3DPoints.row(next_p));

	int edge = NULL;
	pair<int, int> strokeEdge(prev_p, next_p);

	int proj_faceID = -1;
	Eigen::Vector3f bc;
	while(true) {
		if(igl::unproject_onto_mesh(stroke2DPoints.row(next_p), modelview, origin_stroke.viewer.core.proj, origin_stroke.viewer.core.viewport, origin_stroke.get_V(), origin_stroke.get_F(), proj_faceID, bc)) {
			if(proj_faceID == faceID) {
				forward = true;
				return faceID;
			}
		}

		edge = find_next_edge(strokeEdge, edge, faceID, modelview);
		if(edge == NULL) {
			return NULL;
		}

		Eigen::Vector3d v = cutPlane.cross_point(stroke3DPoints.row(EV(edge, 0)), stroke3DPoints.row(EV(edge, 1)));
		PathElement newElement(edge, PathElement::EDGE, v);
		path.push_back(newElement);

		faceID = (EF(edge, 0) == faceID) ? EF(edge, 1) : EF(edge, 0); //get the polygon on the other side of the edge

		//This means that the current stroke point is projected into both "its own" polygon and into the polygon of the next point, WHILE the next point's polygon is across an edge. MUST MEAN that the next point's polygon is on the backside
		if(igl::unproject_onto_mesh(stroke2DPoints.row(prev_p), modelview, origin_stroke.viewer.core.proj, origin_stroke.viewer.core.viewport, origin_stroke.get_V(), origin_stroke.get_F(), proj_faceID, bc)) {
			if(proj_faceID == faceID) {
				forward = false;
				return faceID;
			}
		}
	}
}



int SurfacePath::find_next_edge(pair<int, int> strokeEdge, int prev_edge, int polygon, Eigen::Matrix4f modelview) {
	Eigen::RowVector3d start, end, stroke_start, stroke_end;
	Eigen::MatrixX2d stroke2DPoints = origin_stroke.get_stroke2DPoints();
	Eigen::Vector3d tmp;
	for(int i = 0; i < 3; i++) {
		int edge = FE(polygon, i);
		if(edge != prev_edge) {
			tmp = origin_stroke.get_V().row(EV(edge, 0)).transpose();
			igl::project(tmp, modelview, origin_stroke.viewer.core.proj, origin_stroke.viewer.core.viewport, start);

			tmp = origin_stroke.get_V().row(EV(edge, 1)).transpose();
			igl::project(tmp, modelview, origin_stroke.viewer.core.proj, origin_stroke.viewer.core.viewport, end);

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
bool SurfacePath::edges2D_cross(pair<Eigen::Vector2d, Eigen::Vector2d> edge1, pair<Eigen::Vector2d, Eigen::Vector2d> edge2) {
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

bool SurfacePath::front_facing(int faceID) {
	if((origin_stroke.get_V().row(origin_stroke.get_F()(faceID, 0)) - origin_stroke.viewer.core.camera_center).dot(N_faces.row(faceID)) >= 0) {
		return false;
	} else {
		return true;
	}
}