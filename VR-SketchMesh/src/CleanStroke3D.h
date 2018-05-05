#ifndef _VR_CleanStroke3D_H_
#define _VR_CleanStroke3D_H_
#include <SurfacePath.h>
#include <Eigen/Core>

class CleanStroke3D {
public:
   
private:
	
	std::vector<PathElement> resample_by_length_with_fixes(std::vector<PathElement> path_vertices, double unit_length);

	int find_next_fixed_vertex(std::vector<PathElement> path_vertices, int idx);

	std::vector<PathElement> resample_by_length_sub(std::vector<PathElement> path_vertices, int idx0, int idx1, double unit_length);

	double get_stroke_length(std::vector<PathElement> path_vertices, int start_index, int end_index);

	double vertex_distance(PathElement prev, PathElement next);

};
#endif
