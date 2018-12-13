#ifndef __POLYLINE_PD__
#define __POLYLINE_PD__

#include "mesh.h"
#include "scene.h"

class RotablePolyLine : public Mesh {
public:
	/// initialization
	RotablePolyLine(const RodParameters& p, const VectorX& variables, const int& N);
	static VectorX generatePolyline(const RodParameters& p);
	void buildMasMatrix(const VectorX& variables, const RodParameters& p);

	void Render() override;
	void SetRadius(Scalar r) { _radius = r; }

private:
	void render_lines_polyline();
	void render_orientation_frames();
	void render_points_polyline();
	Scalar _radius;
};

#endif //__POLYLINE_PD__