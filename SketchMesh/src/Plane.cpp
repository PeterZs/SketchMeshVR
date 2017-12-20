#include "Plane.h"
#include <Eigen/Core>

Plane::Plane(Eigen::Vector3d v0, Eigen::Vector3d v1, Eigen::Vector3d v2) {
	Eigen::Vector3d vec0 = v1 - v0;
	Eigen::Vector3d vec1 = v2 - v0;
	normal = vec0.cross(vec1);
	base = v0;
	normal.normalize();
}

Eigen::RowVector3d Plane::cross_point(Eigen::RowVector3d edge_start, Eigen::RowVector3d edge_end) {
	double start_to_surface = signed_distance(edge_start);
	double end_to_surface = signed_distance(edge_end);
	double start_to_end = start_to_surface - end_to_surface;
	Eigen::RowVector3d result = edge_start + (edge_end - edge_start)*(start_to_surface / start_to_end);
	return result;
}

double Plane::signed_distance(Eigen::RowVector3d v) {
	Eigen::Vector3d vec = v.transpose() - base;
	return normal.dot(vec);
}