#ifndef _MESH_H_
#define _MESH_H_
#include <Eigen/Core>
#include <vector>

class Patch;

class Mesh {

public:

    Mesh(Eigen::MatrixXd &V, Eigen::MatrixXi &F, Eigen::VectorXi &vertex_boundary_markers, Eigen::VectorXi &part_of_original_stroke, Eigen::VectorXi &new_mapped_indices, Eigen::VectorXi & sharp_edge, int ID);
	~Mesh();
	Eigen::MatrixXd &V;
	Eigen::MatrixXi &F;
	Eigen::VectorXi &vertex_boundary_markers;
    Eigen::VectorXi &part_of_original_stroke;
	Eigen::VectorXi &new_mapped_indices;
	Eigen::VectorXi &sharp_edge;
	std::vector<Patch*> patches;
	const int ID;
};

#endif
