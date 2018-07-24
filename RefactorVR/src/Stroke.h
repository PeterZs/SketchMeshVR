#ifndef _VR_Stroke_H_
#define _VR_Stroke_H_
#include <igl/opengl/glfw/Viewer.h>
#include <Eigen/Core>

#include <vector>
#include <unordered_map>
class Stroke {
public:
	
	Stroke(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,  igl::opengl::glfw::Viewer &v, int stroke_ID_);
	Stroke(const Stroke& origin);
    Stroke& operator=(Stroke);
	void swap(Stroke & tmp);
	~Stroke();
    
	bool addSegment(Eigen::Vector3f& pos);
	void addSegmentAdd(Eigen::Vector3f & pos);
	void addSegmentCut(Eigen::Vector3f & pos);
	void addSegmentExtrusionBase(Eigen::Vector3f & pos);
	void addSegmentExtrusionSilhouette(Eigen::Vector3f & pos);
	void prepend_first_point();
	void append_final_point();
    void counter_clockwise();
    void strokeReset();
    void update_Positions(Eigen::MatrixXd V);
	bool update_vert_bindings(Eigen::VectorXi & new_mapped_indices, Eigen::VectorXi & edge_boundary_markers, Eigen::VectorXi & sharp_edge, Eigen::VectorXi & vertex_is_fixed);
	void undo_stroke_add(Eigen::VectorXi & edge_boundary_markers, Eigen::VectorXi & sharp_edge, Eigen::VectorXi & vertex_is_fixed);
    bool empty2D() const { return stroke2DPoints.isZero(); }
	bool toLoop();
	std::unordered_map<int, int> generate3DMeshFromStroke(Eigen::VectorXi & edge_boundary_markers, Eigen::VectorXi & vertex_is_fixed, Eigen::MatrixXd & mesh_V, Eigen::MatrixXi & mesh_F);
	bool has_self_intersection();
	bool line_segments_intersect(Eigen::RowVector2d& p1, Eigen::RowVector2d& p2, Eigen::RowVector2d& p3, Eigen::RowVector2d& p4);
	int selectClosestVertex(Eigen::Vector3f pos, double & closest_distance);
	double compute_stroke_diag();
    
    int get_vertex_idx_for_point(int i);
    int get_ID();
    Eigen::MatrixXd get_V() const;
    Eigen::MatrixXi get_F() const;
    Eigen::MatrixXd get_stroke2DPoints() const;
    Eigen::MatrixX3d get3DPoints();
	Eigen::MatrixX3d get3DPointsBack();
	Eigen::MatrixX3d get_hand_pos();
	Eigen::MatrixXi get_hit_faces();
	::std::vector<int> get_closest_vert_bindings();
    void set3DPoints(Eigen::MatrixX3d new_3DPoints);
    void set_closest_vert_bindings(::std::vector<int> new_vert_bindings);


	bool is_loop;
	bool has_points_on_mesh;
	bool has_been_outside_mesh;
	bool has_been_reversed;
	bool starts_on_mesh;
	bool ends_on_mesh;
	Eigen::RowVector3d stroke_color;
	igl::opengl::glfw::Viewer &viewer;
	

	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
	const Eigen::MatrixXd &V;
	const Eigen::MatrixXi &F;
	int stroke_ID; //Non-const for the sake of copy assignment operator 

	Eigen::MatrixX3d stroke3DPoints;
	Eigen::MatrixX3d stroke3DPointsBack;
	Eigen::MatrixXd stroke2DPoints; //Used for early checking if point is new (in screen coordinates)
	Eigen::MatrixXi faces_hit;
	Eigen::MatrixX3d hand_pos_at_draw; //Only used for extrusion base strokes.
	Eigen::RowVector3d cut_stroke_final_point; //Only used for cutting strokes. First point outside of the mesh
	Eigen::RowVectorXd cut_stroke_final_point_2D;
	Eigen::Vector3d pos_before_cut;
	Eigen::Vector3d dir_before_cut;
	Eigen::Vector3d pos_after_cut;
	Eigen::Vector3d dir_after_cut;
	bool prev_point_was_on_mesh; //Used for cutting strokes only. Indicates whether this is the first point outside of the mesh after we've been drawing on the mesh
	Eigen::VectorXd dep;

	std::vector<int> closest_vert_bindings;

	static Eigen::MatrixXd resample_stroke2D(Eigen::MatrixXd & original_stroke2DPoints);
	static void move_to_middle(Eigen::MatrixXd &positions, Eigen::MatrixXd &new_positions);
	static void generate_backfaces(Eigen::MatrixXi &faces, Eigen::MatrixXi &back_faces);
};

#endif
