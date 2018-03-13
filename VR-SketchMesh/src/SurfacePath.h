#ifndef _SurfacePath_H_
#define _SurfacePath_H_
#include <Eigen/Core>
#include "Stroke.h"
#include "Plane.h"

//This class is used for curves that are drawn onto an existing mesh surface, and that need to generate new mesh vertices and edges etc. (whereas the Stroke class simply takes drawn curves as they are)
class PathElement;

class SurfacePath {
public:


	SurfacePath();
	bool create_from_stroke_extrude(const Stroke & stroke);
	bool create_from_stroke_cut(const Stroke& stroke);
	int get_origin_stroke_ID();
    std::vector<PathElement> get_path();
	PathElement& get_path_element(int i);

private:
	int extend_path_cut(int prev_p, int next_p, int faceID, bool & on_front_side, int & edge, bool & first_iter);
	int find_next_edge_cut(std::pair<int, int> strokeEdge, int prev_edge, int polygon, bool on_front_side, ::Plane & cutPlane, Eigen::RowVector3d & start_pos, Eigen::RowVector3d & end_pos, bool first_iter);
	int extend_path_extrude(int prev_p, int next_p, int faceID, Eigen::Matrix4f & modelview);
	int find_next_edge_extrude(int next_p, int prev_p, int prev_edge, int polygon);
  
	Eigen::MatrixX3d create_loop_from_front_and_back(Eigen::MatrixX3d & front_3DPoints, Eigen::MatrixX3d & back_3DPoints);

	std::vector<PathElement> path;
};



class PathElement {
public:
	enum ElementType { EDGE, FACE };
	PathElement(int ID_, ElementType type_, Eigen::Vector3d vertex_) :
		ID(ID_),
		type(type_),
		vertex(vertex_), 
		corresponding_vertex_idx(){
	};

	Eigen::Vector3d get_vertex() { return vertex; };
	ElementType get_type() { return type; };
	int get_ID() { return ID; };
	int get_v_idx() { return corresponding_vertex_idx; };
	void set_v_idx(int idx) { corresponding_vertex_idx = idx; };

private:
	int ID;
	ElementType type;
	Eigen::Vector3d vertex;
	int corresponding_vertex_idx;
};

#endif
