#include <igl/edge_topology.h>
#include <igl/unproject_ray.h>
#include "SurfacePath.h"

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
	Eigen::Matrix4f modelview = stroke.viewervr.get_start_action_view() * stroke.viewervr.corevr.get_model();
	
	int prev_p = stroke.get_stroke2DPoints().rows() - 1;
	int start_p = prev_p;
	int next_p;

	int faceID = origin_stroke->get_hit_faces()(prev_p, 0);
	int start_face = faceID;
	int n = 0;
	Eigen::RowVector3d pt(0, 0, 0);

	int nr3DPoints = origin_stroke->get3DPoints().rows();
	looped_3DPoints = origin_stroke->get3DPoints().topRows(nr3DPoints - 1);//TODO: check how to fix. Need something like we had before in extend path or so (due to order change after counter_clock after to_loop)

	igl::edge_topology(stroke.get_V(), stroke.get_F(), EV, FE, EF);	

	while(true) {
		next_p = n;

		pt = origin_stroke->get3DPoints().row(prev_p);
		//TODO: check if we need to use get_hit_faces like for cut
		PathElement newElement(faceID, PathElement::FACE, pt);
		path.push_back(newElement);

		faceID = extend_path_extrude(prev_p, next_p, faceID, modelview);

		if (next_p == start_p && faceID == start_face) {
			break;
		}

		n = (next_p + 1) % nr3DPoints;
		prev_p = next_p;
	}
	pt = origin_stroke->get3DPoints().row(start_p); //TODO: check why this is needed
	PathElement lastElement(faceID, PathElement::FACE, pt);
	path.push_back(lastElement);
	
}

int SurfacePath::extend_path_extrude(int prev_p, int next_p, int faceID, Eigen::Matrix4f& modelview) {
	Eigen::Vector3d src, dir;
	Eigen::Vector2d tmp = origin_stroke->get_stroke2DPoints().row(prev_p);
	igl::unproject_ray(tmp, modelview, origin_stroke->viewervr.corevr.get_proj(), origin_stroke->viewervr.corevr.viewport, src, dir);
	::Plane cutPlane(src.transpose(), looped_3DPoints.row(prev_p), looped_3DPoints.row(next_p));

	int edge = -1;

	while (true) {

		if (origin_stroke->get_hit_faces()(next_p, 0) == faceID) { //next_p is in same face as prev_p
			return faceID;
		}

		edge = find_next_edge_extrude(next_p, prev_p, edge, faceID);
		if (edge == -1) {
			cout << "This (maybe) shouldn't happen" << endl; //TODO
			return -1;
		}

		//TODO: reuse point that is created by find_next_edge_extrude instead and remove this cutplane here
		Eigen::Vector3d v = cutPlane.cross_point(origin_stroke->get_V().row(EV(edge, 0)), origin_stroke->get_V().row(EV(edge, 1)));
		PathElement newElement(edge, PathElement::EDGE, v);
		path.push_back(newElement);

		faceID = (EF(edge, 0) == faceID) ? EF(edge, 1) : EF(edge, 0); //get the polygon on the other side of the edge
		cout << "new faceID" << faceID << endl;
	}
}

//Find out the edge index of the edge to cross in order to get from polygon to next_p
int SurfacePath::find_next_edge_extrude(int next_p, int prev_p, int prev_edge, int polygon) {
	//polygon is the faceID of prev_p
	int next_faceID = origin_stroke->get_hit_faces()(next_p, 0);

	for (int i = 0; i < 3; i++) {
		int edge = FE(polygon, i);
		if (edge != prev_edge) {
			int other_faceID = (EF(edge, 0) == polygon) ? EF(edge, 1) : EF(edge, 0);
			if (other_faceID == next_faceID) {
				cout << "found adjacent edge" << endl;
				return edge;
			}
		}
	}

	//If we come here, it means that the next strokepoint (next_p) lies in a non-adjacent face (e.g. one or multiple faces away). Only then go for more elaborate checking 
	::Plane cutPlane(origin_stroke->get_hand_pos().row(prev_p), looped_3DPoints.row(prev_p), looped_3DPoints.row(next_p));
	cout << "prev_p " << prev_p << " next_p " << next_p << endl;
	cout << "faces hit: " <<  origin_stroke->get_hit_faces()(prev_p, 0) << "  " << origin_stroke->get_hit_faces()(next_p, 0) << endl;
	cout << "plane points" << endl << origin_stroke->get_hand_pos().row(prev_p) << endl << looped_3DPoints.row(prev_p) << endl << looped_3DPoints.row(next_p) << endl;
	for (int i = 0; i < 3; i++) {
		int edge = FE(polygon, i);
		if (edge != prev_edge) {
			Eigen::RowVector3d cut_point = cutPlane.cross_point(origin_stroke->get_V().row(EV(edge, 0)), origin_stroke->get_V().row(EV(edge, 1)));
			Eigen::RowVector3d seg_vec = looped_3DPoints.row(next_p) - looped_3DPoints.row(prev_p);
			seg_vec.normalize();
			Eigen::RowVector3d its_vec = cut_point - looped_3DPoints.row(prev_p);
			its_vec.normalize();
			if (its_vec.dot(seg_vec) < 0) {
				//Projection points into the opposite direction of line segment
				cout << "test next edge finding for extrude" << endl << seg_vec << endl << its_vec << endl << cut_point << endl;
			}
			else {
				return edge;
			}
		}
	}

	return -1;//shouldn't happen
}

/** Creates a SurfacePath that contains both the original points in stroke, and new points at the locations where stroke segments cross face edges. Also wraps around to the backside of the mesh. Used for cutting **/
void SurfacePath::create_from_stroke_cut(const Stroke & stroke) {
	origin_stroke = new Stroke(stroke);
	path.clear();

	int prev_p = 1; //Start at point 1, because point 0 is defined to be the last point at the beginning of the stroke to lie outside of the mesh
	int start_p = prev_p, next_p;
	int faceID = origin_stroke->get_hit_faces()(prev_p, 0), prev_faceID;

	looped_3DPoints = create_loop_from_front_and_back(origin_stroke->get3DPoints(), origin_stroke->get3DPointsBack());
	int nr_looped_3DPoints = looped_3DPoints.rows();

	int start_face = faceID;
	int n = 2;
	Eigen::RowVector3d pt(0, 0, 0);

	igl::edge_topology(stroke.get_V(), stroke.get_F(), EV, FE, EF);
	bool on_front_side = true, first_iter = true;
	int edge = -1;
	while(true) {
		next_p = n;


		prev_faceID = faceID;
		int prev_p_tmp = (prev_p > nr_looped_3DPoints / 2) ? nr_looped_3DPoints - prev_p : prev_p;
		faceID = origin_stroke->get_hit_faces()(prev_p_tmp, !on_front_side);
		cout << "face and prev face" << faceID << "  " << prev_faceID << endl;
		if (faceID == -1) {
			//We're dealing with one of the 2 off-mesh vertices
			faceID = prev_faceID; //Reset to the previous faceID
		} else {
			pt = looped_3DPoints.row(prev_p);
			PathElement newElement(faceID, PathElement::FACE, pt);
			path.push_back(newElement);
		}

		//bool forward;
		cout << "prev and next" << prev_p << "  " << next_p << endl;
		faceID = extend_path_cut(prev_p, next_p, faceID, on_front_side, edge, first_iter);


		if(next_p == start_p && faceID == start_face) {
			break;
		}

		n = (next_p + 1) % nr_looped_3DPoints;
		prev_p = next_p;
	}

}

//TODO: NOTE that this is now fully adapted to the needs for CUT and not regarding EXTRUDE
/** Determines the moving direction when going from prev_p to next_p and adds new vertices at mesh edges when the segment from prev_p to next_p crosses an edge. **/
int SurfacePath::extend_path_cut(int prev_p, int next_p, int faceID, bool& on_front_side, int& edge, bool& first_iter) {
	//Eigen::MatrixX2d stroke2DPoints = origin_stroke->get_stroke2DPoints();
	//Eigen::MatrixX3d stroke3DPoints = origin_stroke->get3DPoints();
	//int prev_p3D = origin_stroke->has_been_reversed ? prev_p + 1 : prev_p; //These two are needed because when a CW stroke is reversed to a CCW stroke, it generates a discrepancy between the points in the 2DPoints and 3DPoints. 3DPoints is a loop and has vertex 0 also at the end, while 2DPoints doesn't. Fixing this would break many other things
	//int next_p3D = origin_stroke->has_been_reversed ? next_p + 1 : next_p;
	Eigen::RowVector3d plane_point;
	//Don't need to create the cutPlane if prev_p is outside of the mesh, because next_p will be in the same faceID
	int nr_looped_3DPoints = looped_3DPoints.rows();
	int tmp_prev_p = (prev_p > nr_looped_3DPoints / 2) ? nr_looped_3DPoints - prev_p : prev_p;
	if (origin_stroke->get_hit_faces()(tmp_prev_p, 0) == -1) {
	//	plane_point = Eigen::RowVector3d(0, 0, 0);
		int tmp_next_p = nr_looped_3DPoints - next_p;
		plane_point = looped_3DPoints.row(tmp_next_p);
	}
	else {
		int tmp_prev_p = nr_looped_3DPoints - prev_p; 
		plane_point = looped_3DPoints.row(tmp_prev_p); //Get the backside point
	}
	::Plane cutPlane(plane_point, looped_3DPoints.row(prev_p), looped_3DPoints.row(next_p));
	Eigen::RowVector3d start_pos = looped_3DPoints.row(prev_p);
	Eigen::RowVector3d end_pos = looped_3DPoints.row(next_p);

	//int edge = -1;
	cout << "hitfaces" << origin_stroke->get_hit_faces() << endl;

	while(true) {
		next_p = (next_p > nr_looped_3DPoints / 2) ? nr_looped_3DPoints - next_p : next_p; //next_p is used to index into hit_faces, which only has data stored for "the front half" and then has the info for the "back half" in the second column
		prev_p = (prev_p > nr_looped_3DPoints / 2) ? nr_looped_3DPoints - prev_p : prev_p; //prev_p is used to index into hit_faces, which only has data stored for "the front half" and then has the info for the "back half" in the second column
		
		if(origin_stroke->get_hit_faces()(next_p, 0) == faceID || origin_stroke->get_hit_faces()(next_p, 1) == faceID){ //next_p is in same face as prev_p
			return faceID;
		}

		pair<int, int> strokeEdge(prev_p, next_p);
		edge = find_next_edge_cut(strokeEdge, edge, faceID, on_front_side, cutPlane, start_pos, end_pos, first_iter);
		first_iter = false;
		if(edge == -1) {
			cout << "This (maybe) shouldn't happen" << endl; //TODO
			return -1;
		}

		Eigen::Vector3d v = cutPlane.cross_point(origin_stroke->get_V().row(EV(edge, 0)), origin_stroke->get_V().row(EV(edge, 1)));
		PathElement newElement(edge, PathElement::EDGE, v);
		path.push_back(newElement);

		faceID = (EF(edge, 0) == faceID) ? EF(edge, 1) : EF(edge, 0); //get the polygon on the other side of the edge

		if(origin_stroke->get_hit_faces()(prev_p, 0) == faceID || origin_stroke->get_hit_faces()(prev_p, 1) == faceID){// || old_normal.dot(new_normal) < 0){
			cout << " this becomes true" << endl;
            //This means that the current point (prev_p) is projected into both "its own" polygon and into the polygon of next_p, while next_p's polygon is across an edge. Means that next_p's polygon is on the backside
			on_front_side = !on_front_side;
			return faceID;
		}
	} 
}

/** Finds and returns the edge ID of the edge that is being crossed by the segment from strokeEdge.start to strokeEdge.end. Returns -1 if no such edge exists. **/
int SurfacePath::find_next_edge_cut(pair<int, int> strokeEdge, int prev_edge, int polygon, bool on_front_side, ::Plane& cutPlane, Eigen::RowVector3d& start_pos, Eigen::RowVector3d& end_pos, bool first_iter) {
	//polygon is the faceID of prev_p
	int next_faceID;
	if (origin_stroke->get_hit_faces()(strokeEdge.second, !on_front_side) == -1) { //Next point is a point outside of the mesh. Find the edge with prev_p's hit on the other side instead
		int final_prev = (strokeEdge.first > looped_3DPoints.rows() / 2) ? looped_3DPoints.rows() - strokeEdge.first : strokeEdge.first;
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

	
	//If we come here, it means that the next strokepoint (next_p) lies in a non-adjacent face (e.g. one or multiple faces away). Only then go for more elaborate checking 
	cout << polygon << "   " << start_pos << "   " << end_pos << endl;
	double t_val;
	for (int i = 0; i < 3; i++) {
		int edge = FE(polygon, i);
		if (edge != prev_edge) {
			Eigen::RowVector3d cut_point = cutPlane.cross_point(origin_stroke->get_V().row(EV(edge, 0)), origin_stroke->get_V().row(EV(edge, 1)), t_val);
			Eigen::RowVector3d seg_vec = end_pos - start_pos;
			seg_vec.normalize();
			Eigen::RowVector3d its_vec = cut_point - start_pos;
			its_vec.normalize();
			if (t_val>0.99 || t_val<0.01 || (first_iter && its_vec.dot(seg_vec)<0)) { //Force the first extension to go into the right direction
				//Projection points into the opposite direction of line segment
				cout << "test next edge finding for extrude" << endl << seg_vec << endl << its_vec << endl << cut_point << endl;
			}
			else {
				return edge;
			}
		}
	}


	return -1;//shouldn't happen

}

Eigen::MatrixX3d SurfacePath::create_loop_from_front_and_back(Eigen::MatrixX3d& front_3DPoints, Eigen::MatrixX3d& back_3DPoints) {
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
