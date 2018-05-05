#include <CleanStroke3D.h>
#include <vector>

using namespace std;
using namespace igl;

std::vector<PathElement> CleanStroke3D::resample_by_length_with_fixes(std::vector<PathElement> path_vertices, double unit_length) {
	if (path_vertices.size() <= 1) {
		return path_vertices;
	}

	if (path_vertices[0].get_vertex() == path_vertices[path_vertices.size() - 1].get_vertex()) {
		for (int i = 0; i < path_vertices.size(); i++) {
			if (path_vertices[0].fixed) {
				path_vertices.erase(path_vertices.begin() + path_vertices.size() - 1);
				std::vector<PathElement> reordered;
				for (int j = 0; j < path_vertices.size(); j++) {
					reordered.push_back(path_vertices[(i + j) % path_vertices.size()]);
				}
				reordered.push_back(reordered[0]);
				path_vertices = reordered;
				break;
			}
		}
	}

	vector<PathElement> resampled, resample_sub;
	resampled.push_back(path_vertices[0]);

	int idx0 = 0;
	while (true) {
		int idx1 = find_next_fixed_vertex(path_vertices, idx0);
		resample_sub = resample_by_length_sub(path_vertices, idx0, idx1, unit_length);
		resampled.insert(resampled.end(), resample_sub.begin(), resample_sub.end());
		if (idx1 == path_vertices.size() - 1) {
			break;
		}
		idx0 = idx1;
	}
	return resampled;
}

int CleanStroke3D::find_next_fixed_vertex(std::vector<PathElement> path_vertices, int idx) {
	for(int i=idx+1;i<path_vertices.size();i++){
		if (path_vertices[i].fixed) {
			return i;
		}
	}
	return path_vertices.size() - 1;
}

std::vector<PathElement> CleanStroke3D::resample_by_length_sub(std::vector<PathElement> path_vertices, int start_index, int end_index, double unit_length) {
	double length = get_stroke_length(path_vertices, start_index, end_index);
	int n = 1 + (int)(length / unit_length);

	PathElement v0 = path_vertices[start_index];
	PathElement v1 = path_vertices[end_index];

	vector<PathElement> resampled;

	if (length < unit_length) { //Actual stroke is shorter than requested inter-sample length
		resampled.push_back(v1);
		return resampled;
	}

	double total = 0.0, prev_total = 0.0, next_spot = unit_length;
	PathElement prev = v0, next = path_vertices[start_index + 1];
	int index = start_index + 1, count = 0;

	while (true) {
		next = path_vertices[index];
		total += vertex_distance(prev, next);
		while (total >= next_spot) { //The along-path distance to the next path_vertex is bigger than where we would want the next sample, so we create an interpolated sample
			PathElement new_sample = PathElement::interpolate_path_elements(prev, next, (next_spot - prev_total) / (total - prev_total));
			resampled.push_back(new_sample);
			next_spot += unit_length;
			count++;
			if (count == n - 1) {
				break;
			}
		
		}
		if (count == n - 1) {
			break;
		}
		prev = next;
		prev_total = total;
		index++;

		if (index == end_index) {
			break;
		}
	}
	resampled.push_back(v1);
	return resampled;
}

double CleanStroke3D::get_stroke_length(std::vector<PathElement> path_vertices, int start_index, int end_index) {
	double length = 0.0;
	for (int i = start_index; i < end_index; i++) {
		length += vertex_distance(path_vertices[i], path_vertices[i + 1]);
	}
	return length;
}

double CleanStroke3D::vertex_distance(PathElement prev, PathElement next) {
	return (prev.get_vertex() - next.get_vertex()).norm();
}

