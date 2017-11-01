#ifndef IGL_OPENGL2_DRAW_STROKE_H
#define IGL_OPENGL2_DRAW_STROKE_H
#include "../../libigl//include/igl/igl_inline.h"
#include "../../libigl//include/igl/opengl2/gl.h"
#include <vector>
using namespace std;

namespace igl
{
	namespace opengl2
	{
		IGL_INLINE void draw_stroke(const vector<vector<int>> &_stroke);
	}
}

#ifndef IGL_STATIC_LIBRARY
#  include "draw_stroke.cpp"
#endif

#endif
