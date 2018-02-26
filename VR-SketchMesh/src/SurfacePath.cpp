#include <igl/unproject_onto_mesh.h>
#include <igl/edge_topology.h>
#include <igl/project.h>
#include "SurfacePath.h"
#include "Plane.h"

Eigen::MatrixXi EV, FE, EF;
Stroke* origin_stroke;
Eigen::MatrixX3d looped_3DPoints;

using namespace std;

SurfacePath::SurfacePath() {

}

/** Creates a SurfacePath that contains both the original points in stroke, and new points at the locations where stroke segments cross face edges. Won't wrap around to the backside of the mesh (because it will arrive at the first index again before having to switch direction). Used for extrusion **/
void SurfacePath::create_from_stroke_extrude(const Stroke & stroke) {
	origin_stroke = new Stroke(stroke);
	path.clear();
	Eigen::Matrix4f modelview = stroke.viewervr.start_action_view * stroke.viewervr.corevr.model;
	
	//Eigen::Vector3f bc;
	//int faceID = -1;
	int prev_p = stroke.get_stroke2DPoints().rows() - 1;
	int start_p = prev_p;
	int next_p;

	int faceID = origin_stroke->get_hit_faces()(prev_p, 0);
	//igl::unproject_onto_mesh(stroke.get_stroke2DPoints().row(prev_p).cast<float>(), modelview, stroke.viewervr.corevr.proj, stroke.viewervr.corevr.viewport, stroke.get_V(), stroke.get_F(), faceID, bc);
	int start_face = faceID;
	int n = 0;
	Eigen::RowVector3d pt(0, 0, 0);

	igl::edge_topology(stroke.get_V(), stroke.get_F(), EV, FE, EF);

	while(true) {
		next_p = n;

		pt = origin_stroke->get3DPoints().row(prev_p);
		//pt = unproject_onto_polygon(stroke.get_stroke2DPoints().row(prev_p), faceID, modelview);
		PathElement newElement(faceID, PathElement::FACE, pt);
		path.push_back(newElement);

		faceID = extend_path_extrude(prev_p, next_p, faceID, modelview);

		if (next_p == start_p && faceID == start_face) {
			break;
		}

		n = (next_p + 1) % looped_3DPoints.rows();
		prev_p = next_p;
	}
	pt = origin_stroke->get3DPoints().row(start_p); //TODO: check why this is needed
	//pt = unproject_onto_polygon(stroke.get_stroke2DPoints().row(start_p), faceID, modelview);
	PathElement lastElement(faceID, PathElement::FACE, pt);
	path.push_back(lastElement);
	
}

/** Creates a SurfacePath that contains both the original points in stroke, and new points at the locations where stroke segments cross face edges. Also wraps around to the backside of the mesh. Used for cutting **/
void SurfacePath::create_from_stroke_cut(const Stroke & stroke) {
	origin_stroke = new Stroke(stroke);
	path.clear();
	//Eigen::Matrix4f modelview = stroke.viewervr.start_action_view * stroke.viewervr.corevr.model;

	int prev_p = 1; //Start at point 1, because point 0 is defined to be the last point at the beginning of the stroke to lie outside of the mesh
	int start_p = prev_p, next_p;
	int faceID = origin_stroke->get_hit_faces()(prev_p, 0), prev_faceID;

	looped_3DPoints = create_loop_from_front_and_back(origin_stroke->get3DPoints(), origin_stroke->get3DPointsBack());

	int start_face = faceID;
	int n = 2;
	Eigen::RowVector3d pt(0, 0, 0);

	igl::edge_topology(stroke.get_V(), stroke.get_F(), EV, FE, EF);
	bool on_front_side = true;

	while(true) {
		next_p = n;


		pt = looped_3DPoints.row(prev_p);
		prev_faceID = faceID;
		int prev_p_tmp = (prev_p > looped_3DPoints.rows() / 2) ? looped_3DPoints.rows() - prev_p : prev_p;
		faceID = origin_stroke->get_hit_faces()(prev_p_tmp, !on_front_side);
		cout << "check if faceIDs are the same:" << prev_faceID << "  " << faceID << endl;
		if (faceID == -1) {
			//We're dealing with one of the 2 off-mesh vertices
			faceID = prev_faceID; //Reset to the previous faceID
		} else {
			PathElement newElement(faceID, PathElement::FACE, pt);
			path.push_back(newElement);
		}

		//bool forward;
		faceID = extend_path_cut(prev_p, next_p, faceID, on_front_side);


		if(next_p == start_p && faceID == start_face) {
			break;
		}

		n = (next_p + 1) % looped_3DPoints.rows();
		prev_p = next_p;
	}

}

//TODO: NOTE that this is now fully adapted to the needs for CUT and not regarding EXTRUDE
/** Determines the moving direction when going from prev_p to next_p and adds new vertices at mesh edges when the segment from prev_p to next_p crosses an edge. **/
int SurfacePath::extend_path_cut(int prev_p, int next_p, int faceID, bool& on_front_side) {
	//Eigen::MatrixX2d stroke2DPoints = origin_stroke->get_stroke2DPoints();
	//Eigen::MatrixX3d stroke3DPoints = origin_stroke->get3DPoints();
	//int prev_p3D = origin_stroke->has_been_reversed ? prev_p + 1 : prev_p; //These two are needed because when a CW stroke is reversed to a CCW stroke, it generates a discrepancy between the points in the 2DPoints and 3DPoints. 3DPoints is a loop and has vertex 0 also at the end, while 2DPoints doesn't. Fixing this would break many other things
	//int next_p3D = origin_stroke->has_been_reversed ? next_p + 1 : next_p;
	Eigen::RowVector3d plane_point;
	//Don't need to create the cutPlane if prev_p is outside of the mesh, because next_p will be in the same faceID
	int tmp_prev_p = (prev_p > looped_3DPoints.rows() / 2) ? looped_3DPoints.rows() - prev_p : prev_p;
	if (origin_stroke->get_hit_faces()(tmp_prev_p, 0) == -1) {
		plane_point = Eigen::RowVector3d(0, 0, 0);
	}
	else {
		int tmp_prev_p = looped_3DPoints.rows() - prev_p; 
		plane_point = looped_3DPoints.row(tmp_prev_p); //Get the backside point
	}
	::Plane cutPlane(plane_point, looped_3DPoints.row(prev_p), looped_3DPoints.row(next_p));

	int edge = -1;

	while(true) {
		next_p = (next_p > looped_3DPoints.rows() / 2) ? looped_3DPoints.rows() - next_p : next_p; //next_p is used to index into hit_faces, which only has data stored for "the front half" and then has the info for the "back half" in the second column
		prev_p = (prev_p > looped_3DPoints.rows() / 2) ? looped_3DPoints.rows() - prev_p : prev_p; //prev_p is used to index into hit_faces, which only has data stored for "the front half" and then has the info for the "back half" in the second column
		
		if(origin_stroke->get_hit_faces()(next_p, 0) == faceID || origin_stroke->get_hit_faces()(next_p, 1) == faceID){ //next_p is in same face as prev_p
			return faceID;
		}

		pair<int, int> strokeEdge(prev_p, next_p);
		edge = find_next_edge_cut(strokeEdge, edge, faceID, on_front_side);
		if(edge == -1) {
			cout << "This (maybe) shouldn't happen" << endl; //TODO
			return -1;
		}

		Eigen::Vector3d v = cutPlane.cross_point(origin_stroke->get_V().row(EV(edge, 0)), origin_stroke->get_V().row(EV(edge, 1)));
		PathElement newElement(edge, PathElement::EDGE, v);
		path.push_back(newElement);

		faceID = (EF(edge, 0) == faceID) ? EF(edge, 1) : EF(edge, 0); //get the polygon on the other side of the edge
		
		if(origin_stroke->get_hit_faces()(prev_p, 0) == faceID || origin_stroke->get_hit_faces()(prev_p, 1) == faceID){
            //This means that the current point (prev_p) is projected into both "its own" polygon and into the polygon of next_p, while next_p's polygon is across an edge. Means that next_p's polygon is on the backside
			on_front_side = !on_front_side;
			return faceID;
		}
	} 
}

/** Finds and returns the edge ID of the edge that is being crossed by the segment from strokeEdge.start to strokeEdge.end. Returns -1 if no such edge exists. **/
int SurfacePath::find_next_edge_cut(pair<int, int> strokeEdge, int prev_edge, int polygon, bool on_front_side) {
	/*Eigen::MatrixX2d stroke2DPoints = origin_stroke->get_stroke2DPoints();
	Eigen::RowVector2d stroke_start = stroke2DPoints.row(strokeEdge.first).transpose();
	Eigen::RowVector2d stroke_end = stroke2DPoints.row(strokeEdge.second).transpose();
	Eigen::RowVector3d start, end;
	Eigen::RowVector3d tmp;
	for(int i = 0; i < 3; i++) {
		int edge = FE(polygon, i);
		if(edge != prev_edge) {
			tmp = origin_stroke->get_V().row(EV(edge, 0));
			igl::project(tmp, modelview, origin_stroke->viewervr.corevr.proj, origin_stroke->viewervr.corevr.viewport, start);

			tmp = origin_stroke->get_V().row(EV(edge, 1));
			igl::project(tmp, modelview, origin_stroke->viewervr.corevr.proj, origin_stroke->viewervr.corevr.viewport, end);

			if(edges2D_cross({stroke_start, stroke_end}, {start.block(0,0,1,2).transpose(), end.block(0,0,1,2).transpose()})) {
				return edge;
			}
		}
	}
	return -1;*/

	cout << "hit faces: " << endl << origin_stroke->get_hit_faces() << endl << endl;
	//polygon is the faceID of prev_p
	int next_faceID;
	if (origin_stroke->get_hit_faces()(strokeEdge.second, !on_front_side) == -1) { //Next point is a point outside of the mesh. Find the edge with prev_p's hit on the other side instead
		int final_prev = (strokeEdge.first > looped_3DPoints.rows() / 2) ? looped_3DPoints.rows() - strokeEdge.first : strokeEdge.first;
		cout << "stroke edge first" << strokeEdge.first << " " << "final prev" << final_prev << endl;
		next_faceID = origin_stroke->get_hit_faces()(final_prev, on_front_side);
	}
	else {
		int final_next = (strokeEdge.second > looped_3DPoints.rows() / 2) ? looped_3DPoints.rows() - strokeEdge.second : strokeEdge.second; //We don't have hit_face info stored for the backside points (it's in the second index corresponding to the frontside point), so go the frontside
		next_faceID = origin_stroke->get_hit_faces()(final_next, !on_front_side);
	}
	for (int i = 0; i < 3; i++) {
		int edge = FE(polygon, i);
		if (edge != prev_edge) {
			int other_faceID = (EF(edge, 0) == polygon) ? EF(edge, 1) : EF(edge, 0);
			if (other_faceID == next_faceID) {
				return edge;
			}
		}
	}
	return -1;//shouldn't happen

}

int SurfacePath::extend_path_extrude(int prev_p, int next_p, int faceID, Eigen::Matrix4f& modelview) {
	Eigen::Vector3d src, dir;
	Eigen::Vector2d tmp = origin_stroke->get_stroke2DPoints().row(prev_p);
	igl::unproject_ray(tmp, modelview, origin_stroke->viewervr.corevr.proj, origin_stroke->viewervr.corevr.viewport, src, dir);
	::Plane cutPlane(src.transpose(), looped_3DPoints.row(prev_p), looped_3DPoints.row(next_p));

	int edge = -1;
	

	while (true) {
	
		if (origin_stroke->get_hit_faces()(next_p, 0) == faceID) { //next_p is in same face as prev_p
			return faceID;
		}

		edge = find_next_edge_extrude(next_p, edge, faceID);
		if (edge == -1) {
			cout << "This (maybe) shouldn't happen" << endl; //TODO
			return -1;
		}

		Eigen::Vector3d v = cutPlane.cross_point(origin_stroke->get_V().row(EV(edge, 0)), origin_stroke->get_V().row(EV(edge, 1)));
		PathElement newElement(edge, PathElement::EDGE, v);
		path.push_back(newElement);

		faceID = (EF(edge, 0) == faceID) ? EF(edge, 1) : EF(edge, 0); //get the polygon on the other side of the edge
	}
}

//Find out the edge index of the edge to cross in order to get from polygon to next_p
int SurfacePath::find_next_edge_extrude(int next_p, int prev_edge, int polygon) {
	
	//polygon is the faceID of prev_p
	int next_faceID = origin_stroke->get_hit_faces()(next_p, 0);
	
	for (int i = 0; i < 3; i++) {
		int edge = FE(polygon, i);
		if (edge != prev_edge) {
			int other_faceID = (EF(edge, 0) == polygon) ? EF(edge, 1) : EF(edge, 0);
			if (other_faceID == next_faceID) {
				return edge;
			}
		}
	}
	return -1;//shouldn't happen
}

/** Determines whether a pair of 2D segments crosses eachother. Follows principle from https://stackoverflow.com/questions/14176776/find-out-if-2-lines-intersect but slightly different **/
bool SurfacePath::edges2D_cross(pair<Eigen::Vector2d, Eigen::Vector2d> edge1, pair<Eigen::Vector2d, Eigen::Vector2d> edge2) {
	double a0, b0, c0, a1, b1, c1;
	a0 = edge1.first[1] - edge1.second[1];
    b0 = edge1.second[0] - edge1.first[0];
	c0 = edge1.second[1] * edge1.first[0] - edge1.second[0] * edge1.first[1];
	a1 = edge2.first[1] - edge2.second[1];
	b1 = edge2.second[0] - edge2.first[0];
	c1 = edge2.second[1] * edge2.first[0] - edge2.second[0] * edge2.first[1];

	if(((a0*edge2.first[0] + b0*edge2.first[1] + c0)*(a0*edge2.second[0] + b0*edge2.second[1] + c0) <= 0) &&
		((a1*edge1.first[0] + b1*edge1.first[1] + c1)*(a1*edge1.second[0] + b1*edge1.second[1] + c1) <= 0)) {
		return true;
	} else {
		return false;
	}
}

/** Returns whether a face is front facing in the current viewer. **/
bool SurfacePath::front_facing(int faceID) {
	return is_counter_clockwise(faceID);
}

/** Determines whether the vertices of an edge are in counter-clockwise order as seen from the current viewpoint. **/
bool SurfacePath::is_counter_clockwise(int faceID) {
	//TODO: Not sure if we should zero-mean the face vertices first
	Eigen::Matrix4f modelview = origin_stroke->viewervr.start_action_view * origin_stroke->viewervr.corevr.model;
	double total_area = 0.0;
	Eigen::RowVector3d prev, next, pt, center(0, 0, 0);
	Eigen::Vector2d prev2D, next2D;
	for(int i = 0; i < 3; i++) {
		center += origin_stroke->get_V().row(origin_stroke->get_F()(faceID, i));
	}
	center /= 3;

	prev = origin_stroke->get_V().row(origin_stroke->get_F()(faceID, 2)) - center;
	igl::project(prev, modelview, origin_stroke->viewervr.corevr.proj, origin_stroke->viewervr.corevr.viewport, pt);
	prev2D = pt.leftCols(2).transpose();
	for(int i = 0; i < 3; i++) {
		next = origin_stroke->get_V().row(origin_stroke->get_F()(faceID, i)) - center;
		igl::project(next, modelview, origin_stroke->viewervr.corevr.proj, origin_stroke->viewervr.corevr.viewport, pt);
		next2D = pt.leftCols(2).transpose();
		total_area += (prev2D[1] + next2D[1]) * (next2D[0] - prev2D[0]);
		prev2D = next2D;
	}

	if(total_area > 0) {
		return false;
	}
	return true;

}

/** Checks if point v is projected in face as seen from the current viewer. Checks if the point is on the correct side of all face edges. **/
bool SurfacePath::is_projected_inside(Eigen::RowVector2d v, int face, Eigen::Matrix4f modelview) {
	int sign = -1;
	if(!front_facing(face)) {
		sign = 1;
	}

	Eigen::RowVector3d start, end, tmp;
	Eigen::Vector2d vec0, vec1;
	for(int i = 0; i < 3; i++) {
		tmp = origin_stroke->get_V().row((origin_stroke->get_F()(face, i)));
		igl::project(tmp, modelview, origin_stroke->viewervr.corevr.proj, origin_stroke->viewervr.corevr.viewport, start);
		tmp = origin_stroke->get_V().row((origin_stroke->get_F()(face, (i + 1) % 3)));
		igl::project(tmp, modelview, origin_stroke->viewervr.corevr.proj, origin_stroke->viewervr.corevr.viewport, end);

		vec0 = (v - start.block(0, 0, 1, 2)).transpose();
		vec1 = (end.block(0, 0, 1, 2) - start.block(0, 0, 1, 2)).transpose();

		if(cross_prod2D(vec0, vec1)*sign < 0) {
			return false;
		}

	}
	return true;
}

int SurfacePath::cross_prod2D(Eigen::Vector2d vec0, Eigen::Vector2d vec1) {
    return vec0[0] * vec1[1] - vec0[1] * vec1[0];
}

/** Takes a 2D point and unprojects this onto the face with ID faceID. This does not always have to be the first face in the line of sight to the 3D point (as compared to igl::unproject_onto_mesh, which unprojects to the first triangle that is hit). **/
Eigen::Vector3d SurfacePath::unproject_onto_polygon(Eigen::Vector2d point, int faceID, Eigen::Matrix4f modelview) {
	Eigen::Vector3d source, dir;
	Eigen::Vector3d tmp1 = origin_stroke->get_V().row(origin_stroke->get_F()(faceID, 0));
	Eigen::Vector3d tmp2 = origin_stroke->get_V().row(origin_stroke->get_F()(faceID, 1));
	Eigen::Vector3d tmp3 = origin_stroke->get_V().row(origin_stroke->get_F()(faceID, 2));

	double t, u, v;
	igl::unproject_ray(point, modelview, origin_stroke->viewervr.corevr.proj, origin_stroke->viewervr.corevr.viewport, source, dir);
	intersect_triangle1(source.data(), dir.data(), tmp1.data(), tmp2.data(), tmp3.data(), &t, &u, &v);

	return source + t*dir;
}

Eigen::MatrixX3d SurfacePath::create_loop_from_front_and_back(Eigen::MatrixX3d& front_3DPoints, Eigen::MatrixX3d& back_3DPoints) {
	cout << "front " << front_3DPoints << endl << endl << "back " <<back_3DPoints << endl << endl;
	Eigen::MatrixX3d result(front_3DPoints.rows() -1 + back_3DPoints.rows(), 3); //Take all of 3DPointsBack, because we never add back-points for the points that are off the mesh
	result << front_3DPoints.topRows(front_3DPoints.rows() - 1), back_3DPoints.colwise().reverse().eval();
	return result;
}

vector<PathElement> SurfacePath::get_path() {
	return path;
}

int SurfacePath::get_origin_stroke_ID() {
	return origin_stroke->get_ID();
}

PathElement & SurfacePath::get_path_element(int i) {
	return path[i];
}
