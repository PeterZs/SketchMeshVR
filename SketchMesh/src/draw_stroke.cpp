#include "draw_stroke.h"
#include "../../libigl//include/igl/opengl2/gl.h"
#include <vector>

using namespace std;
IGL_INLINE void igl::opengl2::draw_stroke(const vector<vector<int>> & _stroke) {
	using namespace std;
	glBegin(GL_LINE_STRIP);
	for (int i = 0; i < _stroke.size(); i++) {
		glVertex2f(_stroke[i][0], _stroke[i][1]);
	}
	glEnd();
}