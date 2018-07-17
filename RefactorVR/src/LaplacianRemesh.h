#ifndef _Laplacian_Remesh_H_
#define _Laplacian_Remesh_H_
#include "Mesh.h"
#include "SurfacePath.h"

class LaplacianRemesh {

public:
	static Eigen::VectorXi remesh_cut_remove_inside(Mesh & m, SurfacePath & surface_path, Eigen::Matrix4f model, Eigen::Matrix4f view, Eigen::Matrix4f proj, Eigen::Vector4f viewport, bool & remesh_success, int cut_clicked_face);
	static Eigen::VectorXi remesh_extrusion_remove_inside(Mesh & m, SurfacePath & surface_path, Eigen::Matrix4f model, Eigen::Matrix4f view, Eigen::Matrix4f proj, Eigen::Vector4f viewport, bool& remesh_success);
	static bool remesh_open_path(Mesh & m, Stroke & open_path_stroke);

	static bool remesh_cutting_path(Mesh & m, Stroke & cut_path_stroke);
	static Eigen::VectorXi remesh(Mesh & m, SurfacePath & surface_path, Eigen::Matrix4f model, Eigen::Matrix4f view, Eigen::Matrix4f proj, Eigen::Vector4f viewport, bool & remesh_success, int cut_clicked_face = -1);

private:
	static bool is_front_loop;
	static bool remove_inside_faces;

	static Eigen::MatrixXi EV, FE, EF;
	static std::vector<std::vector<int>> VV;

	static int find_closest(std::vector<int> vertices, Eigen::Vector3d base, Mesh & m);
	static void stitch(std::vector<int> path_vertices, std::vector<int> boundary_vertices, Mesh & m, bool reverse);
	static bool is_counter_clockwise_boundaries(Eigen::MatrixXd boundary_points, Eigen::Matrix4f modelview, Eigen::Matrix4f proj, Eigen::Vector4f viewport, Eigen::RowVector3d mean_viewpoint, bool cut);
	static double compute_average_distance_between_onPolygon_vertices(std::vector<PathElement> path);
	static double compute_average_length_of_crossing_edges(std::vector<PathElement> path, Eigen::MatrixXd startV, Eigen::MatrixXi startEV);
	static void update_add_path_points_and_bindings(Stroke & stroke, SurfacePath & surface_path);
	static std::vector<int> sort_boundary_vertices(Eigen::Vector3d start_vertex, std::vector<int> boundary_vertices, Mesh & m);
	static std::vector<int> reorder(std::vector<int> boundary_vertices, Eigen::Vector3d start_v, Mesh & m);
	static Eigen::RowVector3d compute_mean_viewpoint(Mesh & m, std::vector<int> inner_boundary_vertices);
	static void propagate_dirty_faces(int face, std::vector<bool>& dirty_face);
	static void update_sharp_edges(Mesh & m, Eigen::MatrixXi & all_sharpEV);
	static void update_mesh_values(Mesh & m, Eigen::MatrixXd path, int stroke_ID, int new_mapped_start, bool hold_back_due_to_loop);
	static void update_face_indices(Mesh & m);



};

#endif