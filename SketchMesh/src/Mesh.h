#ifndef _MESH_H_
#define _MESH_H_
#include <Eigen/Core>
#include <Eigen/Geometry>

class Mesh {

public:

    Mesh(Eigen::MatrixXd &V, Eigen::MatrixXi &F, Eigen::VectorXi &vertex_boundary_markers, Eigen::VectorXi &part_of_original_stroke, int ID);
	~Mesh();
	Eigen::MatrixXd &V;
	Eigen::MatrixXi &F;
	Eigen::VectorXi &vertex_boundary_markers;
    Eigen::VectorXi &part_of_original_stroke;
	const int ID;
};

#endif
