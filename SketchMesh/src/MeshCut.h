#ifndef _Mesh_Cut_H_
#define _Mesh_Cut_H_
#include <Eigen/Core>
#include "Stroke.h"
#include "Mesh.h"
#include "SurfacePath.h"

class MeshCut {

public:
	static void cut(Eigen::MatrixXd & V, Eigen::MatrixXi & F, Eigen::VectorXi & vertex_boundary_markers, Eigen::VectorXi & part_of_original_stroke, Stroke & stroke);


private:
	static int prev_vertex_count;
	static int ID;
	static Eigen::MatrixXi EV, FE, EF;

	static void cut_main(Mesh& m, SurfacePath& surface_path, Stroke& stroke);

	//static Eigen::VectorXi remesh_cut_remove_inside(Mesh & m, SurfacePath & surface_path);
	//static void propagate_dirty_faces(int face, std::vector<bool>& dirty_face);

	//static std::vector<int> sort_boundary_vertices(Eigen::Vector3d start_vertex, std::vector<int> boundary_vertices, Mesh& m);
	//static void stitch(std::vector<int> path_vertices, std::vector<int> outer_boundary_vertices, Mesh& m);
	//static std::vector<int> MeshCut::reorder(std::vector<int> boundary_vertices, Eigen::Vector3d start_v, Mesh& m);
	//static int find_closest(std::vector<int> vertices, Eigen::Vector3d base, Mesh& m);
	//static void reverse_path(Eigen::MatrixXd path_vertices);
	//static void reverse_path(std::vector<int> path_vertices);
	//static double compute_average_length_of_crossing_edges(std::vector<PathElement> path, Mesh& m);
	//static double compute_average_distance_between_onPolygon_vertices(std::vector<PathElement> path, Mesh & m);
	//static double get_length(Eigen::MatrixX3d path_vertices);
	static Eigen::MatrixXd resample_by_length_with_fixes(std::vector<int> path_vertices, Mesh & m, double unit_length);
	//static Eigen::MatrixX3d resample_stroke(Eigen::MatrixX3d & original_stroke3DPoints);
	//static void move_to_middle(Eigen::MatrixX3d & positions, Eigen::MatrixX3d & new_positions);
	static void mesh_open_hole(Eigen::VectorXi & boundary_vertices, Mesh & m);


	static void post_cut_update_points(Stroke& stroke, SurfacePath& surface_path);

};


#endif

