#ifndef _Laplacian_Remesh_H_
#define _Laplacian_Remesh_H_
#include <Eigen/Core>
#include "Mesh.h"
#include "SurfacePath.h"

class LaplacianRemesh {

public:
	static Eigen::VectorXi remesh_cut_remove_inside(Mesh & m, SurfacePath & surface_path, Eigen::Matrix4f model, Eigen::Matrix4f view, Eigen::Matrix4f proj, Eigen::Vector4f viewport);
	static Eigen::VectorXi remesh_extrusion_remove_inside(Mesh & m, SurfacePath & surface_path, Eigen::Matrix4f model, Eigen::Matrix4f view, Eigen::Matrix4f proj, Eigen::Vector4f viewport);
	static Eigen::VectorXi remesh(Mesh & m, SurfacePath & surface_path, Eigen::Matrix4f model, Eigen::Matrix4f view, Eigen::Matrix4f proj, Eigen::Vector4f viewport);

private:
	static bool is_front_loop;
	static bool remove_inside_faces;

	static Eigen::MatrixXi EV, FE, EF;
	static std::vector<std::vector<int>> VV;


	static void propagate_dirty_faces(int face, std::vector<bool>& dirty_face);

	static std::vector<int> sort_boundary_vertices(Eigen::Vector3d start_vertex, std::vector<int> boundary_vertices, Mesh & m);

	static int find_closest(std::vector<int> vertices, Eigen::Vector3d base, Mesh & m);

	static void stitch(std::vector<int> path_vertices, std::vector<int> boundary_vertices, Mesh & m);

	static std::vector<int> reorder(std::vector<int> boundary_vertices, Eigen::Vector3d start_v, Mesh & m);

	static void reverse_path(std::vector<int> path_vertices);

	static double compute_average_length_of_crossing_edges(std::vector<PathElement> path, Mesh & m);

	static double compute_average_distance_between_onPolygon_vertices(std::vector<PathElement> path, Mesh & m);

	static double get_length(Eigen::MatrixX3d path_vertices);

	static Eigen::MatrixX3d resample_stroke(Eigen::MatrixX3d & original_stroke3DPoints);

	static void move_to_middle(Eigen::MatrixX3d & positions, Eigen::MatrixX3d & new_positions);

	static bool is_counter_clockwise_boundaries(Eigen::MatrixXd boundary_points, Eigen::Matrix4f modelview, Eigen::Matrix4f proj, Eigen::Vector4f viewport, Eigen::RowVector3d mean_viewpoint, bool cut);

};

#endif