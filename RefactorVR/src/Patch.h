#ifndef _PATCH_H_
#define _PATCH_H_

#include "Mesh.h"


class Patch {

public:
	Patch();
	~Patch();

	static std::vector<Patch*> init_patches(Mesh& h);
	static void propagate_patch(Patch& patch, int face, Eigen::VectorXi& faces);

	Mesh mesh;

private:


};
#endif
