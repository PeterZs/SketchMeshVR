#ifndef _Plane_H_
#define _Plane_H_
#include <Eigen/Core>



class Plane {
public:
	Eigen::Vector3d base;
	Eigen::Vector3d normal;

	Plane(Eigen::Vector3d v0, Eigen::Vector3d v1, Eigen::Vector3d v2);
	Eigen::RowVector3d cross_point(Eigen::RowVector3d edge_start, Eigen::RowVector3d edge_end);

private:
	double signed_distance(Eigen::RowVector3d v);

};
#endif