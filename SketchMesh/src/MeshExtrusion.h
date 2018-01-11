#ifndef _Mesh_Extrusion_H_
#define _Mesh_Extrusion_H_
#include <Eigen/Core>
#include "Stroke.h"
#include "Mesh.h"
#include "SurfacePath.h"

class MeshExtrusion {

public:
	static void extrude(Eigen::MatrixXd &V, Eigen::MatrixXi &F, Eigen::VectorXi &vertex_boundary_markers, Eigen::VectorXi &part_of_original_stroke, Stroke& base, Stroke& silhouette);


private:
	static int prev_vertex_count;
	static int ID;

	static void extrude_main(Mesh & m, SurfacePath & surface_path, Stroke & stroke);

};

#endif