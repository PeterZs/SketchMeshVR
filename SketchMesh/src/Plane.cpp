#include "Plane.h"
#include <Eigen/Core>
#include <Eigen/Geometry>

Plane::Plane(Eigen::Vector3d v0, Eigen::Vector3d v1, Eigen::Vector3d v2) {
	Eigen::Vector3d vec0 = v1 - v0;
	Eigen::Vector3d vec1 = v2 - v0;
	normal = vec0.cross(vec1);
	base = v0;
	normal.normalize();
}

Plane::Plane(Eigen::Vector3d point, Eigen::Vector3d normal_):
base(point),
normal(normal_) {
} 

Eigen::RowVector3d Plane::cross_point(Eigen::RowVector3d edge_start, Eigen::RowVector3d edge_end) {
	double start_to_surface = signed_distance(edge_start);
	double end_to_surface = signed_distance(edge_end);
	double start_to_end = start_to_surface - end_to_surface;
	Eigen::RowVector3d result = edge_start + (edge_end - edge_start)*(start_to_surface / start_to_end);
	return result;
}

Eigen::RowVector3d Plane::cross_point_ray(Eigen::Vector3d ray_source, Eigen::Vector3d ray_dir) {
	double plane_to_source = signed_distance(ray_source);
	double cos_ = cos(-normal, ray_dir);
	return ray_source + ray_dir*(plane_to_source/cos_);
}

double Plane::signed_distance(Eigen::RowVector3d v) {
	Eigen::Vector3d vec = v.transpose() - base;
	return normal.dot(vec);
}

//TODO: convince myself that this works. comes from teddy. rename and move to better place? Comes from teddy.vector3
double Plane::cos(Eigen::Vector3d u, Eigen::Vector3d v) {
	double length = sqrt((u.dot(u))*(v.dot(v)));
	if(length > 0) {
		return u.dot(v) / length;
	} else {
		return 0;
	}
}