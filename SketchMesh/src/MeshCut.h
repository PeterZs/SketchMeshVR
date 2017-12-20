#ifndef _Mesh_Cut_H_
#define _Mesh_Cut_H_
#include <Eigen/Core>
#include "Stroke.h"

class MeshCut {

public:
	static void cut(Eigen::MatrixXd & V, Eigen::MatrixXi & F, Eigen::VectorXi & vertex_boundary_markers, Eigen::VectorXi & part_of_original_stroke, Stroke & stroke);


private:
	static int prev_vertex_count;
	static int ID;


};


#endif

