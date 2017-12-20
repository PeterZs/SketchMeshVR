#ifndef _SurfacePath_H_
#define _SurfacePath_H_
#include <Eigen/Core>
#include "Stroke.h"

//This class  is used for curves that are drawn onto an existing mesh surface, and that need to generate new mesh vertices and edges etc. (whereas the Stroke class simply takes drawn curves as they are)

class SurfacePath {
public:


	SurfacePath();
	void create_from_stroke(Stroke& stroke);


private:
	int extend_path(int prev_p, int next_p, int faceID, bool forward, bool front_facing_only, Eigen::Matrix4f modelview);
	int find_next_edge(std::pair<int, int> strokeEdge, int prev_edge, int polygon, Eigen::Matrix4f modelview);
	bool edges2D_cross(std::pair<Eigen::Vector2d, Eigen::Vector2d> edge1, std::pair<Eigen::Vector2d, Eigen::Vector2d> edge2);
	bool front_facing(int faceID);
	std::vector<PathElement> path;
};



class PathElement {
public:
	enum ElementType { EDGE, FACE };
	PathElement(int ID_, ElementType type_, Eigen::Vector3d vertex_) :
		ID(ID_),
		type(type_),
		vertex(vertex_){
	};

private:
	int ID;
	ElementType type;
	Eigen::Vector3d vertex;
};

#endif