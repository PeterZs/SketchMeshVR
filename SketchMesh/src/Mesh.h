#ifndef _MESH_
#define _MESH_
#include <Eigen/Core>
#include <Eigen/Geometry>

class Mesh {
	//Contains (a part of) the vertices and faces that are visible on the screen. TODO: this is because what is displayed might consist of several separate "meshes" that are individually optimized
public:

	Mesh(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F);
	~Mesh();
	const Eigen::MatrixXd &V;
	const Eigen::MatrixXi &F;
};

#endif _MESH_