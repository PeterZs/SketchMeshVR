#ifndef _MESH_H_
#define _MESH_H_
#include <Eigen/Core>
#include <Eigen/Geometry>

class Mesh {
	//Contains (a part of) the vertices and faces that are visible on the screen. TODO: this is because what is displayed might consist of several separate "meshes" that are individually optimized
public:

	Mesh(Eigen::MatrixXd &V, Eigen::MatrixXi &F, const Eigen::VectorXi &vertex_boundary_markers, int ID);
	~Mesh();
	Eigen::MatrixXd &V;
	Eigen::MatrixXi &F;
	const Eigen::VectorXi &vertex_boundary_markers;
	const int ID;
};

#endif