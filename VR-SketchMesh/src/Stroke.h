#ifndef _VR_Stroke_H_
#define _VR_Stroke_H_
#include <igl/viewer/VR_Viewer.h>
#include <Eigen/Core>

#include <vector>
#include <unordered_map>
class Stroke {
public:
	
	Stroke(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,  igl::viewer::VR_Viewer &v, int stroke_ID_);
	Stroke(const Stroke& origin);
    Stroke& operator=(Stroke);
	void swap(Stroke & tmp);
	~Stroke();
    
	void strokeAddSegment(Eigen::Vector3f& pos);
	bool strokeAddSegmentAdd(Eigen::Vector3f & pos);
	void strokeAddSegmentCut(Eigen::Vector3f & pos);
	void strokeAddSegmentExtrusionBase(Eigen::Vector3f & pos);
	void strokeAddSegmentExtrusionSilhouette(Eigen::Vector3f & pos);
	void prepend_first_point();
	void append_final_point();
    void counter_clockwise();
    void strokeReset();
    void snap_to_vertices(Eigen::VectorXi & vertex_boundary_markers);
    void undo_stroke_add(Eigen::VectorXi & vertex_boundary_markers);
    void update_Positions(Eigen::MatrixXd V);
    bool update_vert_bindings(Eigen::VectorXi & new_mapped_indices, Eigen::VectorXi & vertex_boundary_markers);
    bool empty2D() const { return stroke2DPoints.isZero(); }
	bool toLoop();
    std::unordered_map<int, int> generate3DMeshFromStroke(Eigen::VectorXi &vertex_boundary_markers, Eigen::VectorXi &part_of_original_stroke);
    int selectClosestVertex(Eigen::Vector3f& pos, double& closest_distance);
    double compute_stroke_diag();
    
    int get_vertex_idx_for_point(int i);
    int get_ID();
    Eigen::MatrixXd get_V() const;
    Eigen::MatrixXi get_F() const;
    Eigen::MatrixXd get_stroke2DPoints() const;
    Eigen::MatrixX3d get3DPoints();
	Eigen::MatrixX3d get3DPointsBack();
	Eigen::MatrixXi get_hit_faces();
	std::vector<int> get_closest_vert_bindings();
    void set3DPoints(Eigen::MatrixX3d new_3DPoints);
    void set_closest_vert_bindings(std::vector<int> new_vert_bindings);


	bool is_loop;
	bool has_points_on_mesh;
	bool has_been_reversed;
	Eigen::RowVector3d stroke_color;
	igl::viewer::VR_Viewer &viewervr;
	

	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
	const Eigen::MatrixXd &V;
	const Eigen::MatrixXi &F;
	int stroke_ID; //Non-const for the sake of copy assignment operator 

	Eigen::MatrixX3d stroke3DPoints;
	Eigen::MatrixX3d stroke3DPointsBack;
	Eigen::MatrixXd stroke2DPoints; //Used for early checking if point is new (in screen coordinates)
	Eigen::MatrixXi faces_hit;
	Eigen::RowVector3d cut_stroke_final_point; //Only used for cutting strokes. First point outside of the mesh
	Eigen::RowVectorXd cut_stroke_final_point_2D;
	Eigen::Vector3d pos_before_cut;
	Eigen::Vector3d dir_before_cut;
	Eigen::Vector3d pos_after_cut;
	Eigen::Vector3d dir_after_cut;
	bool just_came_from_mesh; //Used for cutting strokes only. Indicates whether this is the first point outside of the mesh after we've been drawing on the mesh
	Eigen::VectorXd dep;

	std::vector<int> closest_vert_bindings;

	static Eigen::MatrixXd resample_stroke2D(Eigen::MatrixXd & original_stroke2DPoints);
	static void move_to_middle(Eigen::MatrixXd &positions, Eigen::MatrixXd &new_positions);
	static void generate_backfaces(Eigen::MatrixXi &faces, Eigen::MatrixXi &back_faces);
};

#endif
