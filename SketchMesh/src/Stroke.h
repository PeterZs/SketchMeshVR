#ifndef _Stroke_H_
#define _Stroke_H_
#include <igl/viewer/Viewer.h>
#include <vector>
#include <unordered_map>

class Stroke {
public:
	
	Stroke(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,  igl::viewer::Viewer &v, int stroke_ID_);
	Stroke(const Stroke& origin);
	Stroke& operator=(Stroke);
	void swap(Stroke & tmp);
	~Stroke();
	bool empty2D() const { return stroke2DPoints.isZero(); }
	bool empty3D() const { return stroke3DPoints.isZero(); }
	void strokeAddSegment(int mouse_x, int mouse_y);
	bool strokeAddSegmentAdd(int mouse_x, int mouse_y);
	bool strokeAddSegmentCut(int mouse_x, int mouse_y);
	void strokeAddSegmentExtrusion(int mouse_x, int mouse_y);
	bool toLoop();
	void strokeReset();
    std::unordered_map<int, int> generate3DMeshFromStroke(Eigen::VectorXi &vertex_boundary_markers, Eigen::VectorXi &part_of_original_stroke);
	//double total_stroke_length();
	Eigen::MatrixX3d get3DPoints();
	int get_vertex_idx_for_point(int i);
	int get_ID();
	void prepare_for_cut();
	std::vector<int> get_closest_vert_bindings();
	int selectClosestVertex(int mouse_x, int mouse_y, double& closest_distance);
	double compute_stroke_diag();
	void update_Positions(Eigen::MatrixXd V);
	void snap_to_vertices(Eigen::VectorXi & vertex_boundary_markers);

	void mirror_on_backside(Eigen::VectorXi & vertex_boundary_markers, std::unordered_map<int, int> backside_vertex_map);

	bool is_loop;
	bool has_points_on_mesh;
	Eigen::RowVector3d stroke_color;
	igl::viewer::Viewer &viewer;
	Eigen::MatrixXd get_V();

	Eigen::MatrixXi get_F();

	Eigen::MatrixX2d get_stroke2DPoints();

	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
	const Eigen::MatrixXd &V;
	const Eigen::MatrixXi &F;
	int stroke_ID; //Non-const for the sake of copy assignment operator 

	Eigen::MatrixX3d stroke3DPoints; //Used for screen output
	Eigen::MatrixX2d stroke2DPoints; //Used for early checking if point is new (in screen coordinates)
	Eigen::MatrixXi stroke_edges;
	double dep = -1;

	std::vector<int> closest_vert_bindings;
	std::vector<int> added_stroke_final_vertices;

	void counter_clockwise();
	static Eigen::MatrixX2d resample_stroke(Eigen::MatrixX2d & original_stroke2DPoints);
	static void move_to_middle(Eigen::MatrixX2d &positions, Eigen::MatrixX2d &new_positions);
	static void generate_backfaces(Eigen::MatrixXi &faces, Eigen::MatrixXi &back_faces);

	int extend_path(int prev_p, int next_p, int faceID, Eigen::VectorXi forward, bool front_facing_only, Eigen::Matrix4f modelview);

	int find_next_edge(std::pair<int, int> strokeEdge, int prev_edge, int polygon, Eigen::Matrix4f modelview);

	bool edges2D_cross(std::pair<Eigen::Vector2d, Eigen::Vector2d> edge1, std::pair<Eigen::Vector2d, Eigen::Vector2d> edge2);


};

#endif
