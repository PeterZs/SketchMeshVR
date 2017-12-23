#ifndef _Plane_H_
#define _Plane_H_
#include <Eigen/Core>



class Plane {
public:


	Plane(Eigen::Vector3d v0, Eigen::Vector3d v1, Eigen::Vector3d v2);
	Plane(Eigen::Vector3d point, Eigen::Vector3d normal);
	Eigen::RowVector3d cross_point(Eigen::RowVector3d edge_start, Eigen::RowVector3d edge_end);
	Eigen::RowVector3d cross_point_ray(Eigen::Vector3d ray_source, Eigen::Vector3d ray_dir);

	Eigen::Vector3d base;
	Eigen::Vector3d normal;

private:
	double signed_distance(Eigen::RowVector3d v);

	double cos(Eigen::Vector3d u, Eigen::Vector3d v);

};
#endif