#ifndef _MESH_H_
#define _MESH_H_
#include <Eigen/Core>
#include <Eigen/Geometry>

class Mesh {

public:

    Mesh(Eigen::MatrixXd &V, Eigen::MatrixXi &F, Eigen::VectorXi &vertex_boundary_markers, Eigen::VectorXi &part_of_original_stroke, Eigen::VectorXi &new_mapped_indices, int ID);
	~Mesh();
	Eigen::MatrixXd &V;
	Eigen::MatrixXi &F;
	Eigen::VectorXi &vertex_boundary_markers;
    Eigen::VectorXi &part_of_original_stroke;
	Eigen::VectorXi &new_mapped_indices;
	const int ID;
};

#endif
