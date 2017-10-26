#include <iostream>
#include <igl/readOFF.h>
#include <igl/viewer/Viewer.h>
#include <igl/vertex_triangle_adjacency.h>
#include <igl/adjacency_list.h>
#include <igl/per_face_normals.h>
#include <igl/per_vertex_normals.h>
#include <igl/per_corner_normals.h>
#include <igl/facet_components.h>
#include <igl/jet.h>
#include <igl/barycenter.h>
#include <cmath>
#include <igl/triangle_triangle_adjacency.h>
/*** insert any libigl headers here ***/

using namespace std;
using Viewer = igl::viewer::Viewer;

// Vertex array, #V x3
Eigen::MatrixXd V;
// Face array, #F x3
Eigen::MatrixXi F;
// Per-face normal array, #F x3
Eigen::MatrixXd FN;
// Per-vertex normal array, #V x3
Eigen::MatrixXd VN;
// Per-corner normal array, (3#F) x3
Eigen::MatrixXd CN;
// Vectors of indices for adjacency relations
std::vector<std::vector<int> > VF, VFi, VV;
// Integer vector of component IDs per face, #F x1
Eigen::VectorXi cid;
// Per-face color array, #F x3
Eigen::MatrixXd component_colors_per_face;

bool callback_key_down(Viewer& viewer, unsigned char key, int modifiers) {
	if (key == '1') {
		viewer.data.clear();
		viewer.data.set_mesh(V, F);
		// Add your code for computing vertex to face relations here;
		// store in VF,VFi.
		igl::vertex_triangle_adjacency(V.rows(), F, VF, VFi);
		std::cout << "Vertex" << "      " << "Incident faces" << std::endl;
		for (int i = 0; i<VFi.size(); i++) {
			std::cout << i << "         ";
			std::cout << VF[i][0];
			for (int j = 1; j<VF[i].size(); j++) {
				std::cout << ", " << VF[i][j];
			}
			std::cout << std::endl;
		}

	}

	if (key == '2') {
		viewer.data.clear();
		viewer.data.set_mesh(V, F);
		// Add your code for computing vertex to vertex relations here:
		// store in VV.
		igl::adjacency_list(F, VV);
		std::cout << "Vertex" << "      " << "Adjacent vertices" << std::endl;
		for (int i = 0; i<VV.size(); i++) {
			std::cout << i << "         ";
			std::cout << VV[i][0];
			for (int j = 1; j<VV[i].size(); j++) {
				std::cout << ", " << VV[i][j];
			}
			std::cout << std::endl;
		}

	}

	if (key == '3') {
		viewer.data.clear();
		viewer.data.set_mesh(V, F);
		FN.setZero(F.rows(), 3);
		// Add your code for computing per-face normals here: store in FN.
		igl::per_face_normals(V, F, FN);
		// Set the viewer normals.
		viewer.data.set_normals(FN);
	}

	if (key == '4') {
		viewer.data.clear();
		viewer.data.set_mesh(V, F);
		// Add your code for computing per-vertex normals here: store in VN.
		igl::per_vertex_normals(V, F, VN);
		// Set the viewer normals.
		viewer.data.set_normals(VN);
	}

	if (key == '5') {
		viewer.data.clear();
		viewer.data.set_mesh(V, F);
		// Add your code for computing per-corner normals here: store in CN.
		int threshold = 50;
		igl::per_corner_normals(V, F, threshold, CN);
		//Set the viewer normals
		viewer.data.set_normals(CN);
	}

	if (key == '6') {
		viewer.data.clear();
		viewer.data.set_mesh(V, F);
		component_colors_per_face.setZero(F.rows(), 3);
		// Add your code for computing per-face connected components here:
		// store the component labels in cid.
		igl::facet_components(F, cid);
		// Compute colors for the faces based on components, storing them in
		// component_colors_per_face.
		igl::jet(cid, true, component_colors_per_face);
		// Set the viewer colors
		viewer.data.set_colors(component_colors_per_face);

		std::cout << "Number of connected components: " << cid.maxCoeff() + 1 << std::endl;
		std::cout << "Number of faces per component: " << std::endl;
		for (int i = 0; i<cid.maxCoeff() + 1; i++) {
			std::cout << "Component " << i + 1 << "   " << (cid.array() == i).count() << endl;
		}
		std::cout << endl;
	}

	if (key == '7') {
		Eigen::MatrixXd Vout = V;
		Eigen::MatrixXi Fout = F;
		// Add your code for sqrt(3) subdivision here.
		Eigen::MatrixXd BC;
		igl::barycenter(V, F, BC); //Compute barycenter of each face and store in BC
		Eigen::MatrixXd Vpp(V.rows() + BC.rows(), V.cols()); //Create V''
		Vpp << V, BC; //Put V and BC into V''

		Eigen::MatrixXi Fpp(3 * F.rows(), F.cols());
		for (int i = 0; i<F.rows(); i++) {
			Fpp.row(i * 3 + 0) << (V.rows() + i), F(i, 1), F(i, 2);
			Fpp.row(i * 3 + 1) << F(i, 0), (V.rows() + i), F(i, 2);
			Fpp.row(i * 3 + 2) << F(i, 0), F(i, 1), (V.rows() + i);
		}

		Fout = Fpp;
		Vout = Vpp;

		//Step 2
		Eigen::MatrixXd Vp(Vpp.rows(), Vpp.cols()); //Create V'
		igl::adjacency_list(F, VV);

		for (int i = 0; i<V.rows(); i++) {
			double ni = VV[i].size();
			double an = (4.0 - 2 * cos(2 * M_PI / ni)) / 9.0;
			Eigen::MatrixXd Vsum = Eigen::MatrixXd::Zero(1, 3);
			for (int j = 0; j<ni; j++) {
				Vsum += V.row(VV[i][j]);
			}
			Eigen::MatrixXd pi = (1.0 - an)*V.row(i) + (an / ni) * Vsum;
			Vp.row(i) << pi;
		}
		for (int i = V.rows(); i<V.rows() + BC.rows(); i++) {
			Vp.row(i) << BC.row(i - V.rows());
		}


		//Step 3
		std::vector<std::vector<std::vector<int>>> TT, TTi;
		igl::triangle_triangle_adjacency(Fpp, TT, TTi);

		Eigen::MatrixXi Fp(3 * F.rows(), F.cols());

		//CHECK BOUNDARY CONDITIONS
		for (int i = 0; i<Fpp.rows(); i++) {//Go over the current faces
			if (TT[i][0].size()>0 && TT[i][1].size()>0 && TT[i][2].size()>0) {//Non-boundary face
				if (Fpp(i, 0) >= V.rows() && Fpp(TT[i][0][0], TTi[i][0][0]) >= V.rows()) {//If both vertices are non-original, connect them
					Fp.row(i) << Fpp(i, 0), Fpp(i, 1), Fpp(TT[i][0][0], TTi[i][0][0]);
				}
				if (Fpp(i, 1) >= V.rows() && Fpp(TT[i][1][0], TTi[i][1][0]) >= V.rows()) {//If both vertices are non-original, connect them
					Fp.row(i) << Fpp(i, 1), Fpp(i, 2), Fpp(TT[i][1][0], TTi[i][1][0]);
				}
				if (Fpp(i, 2) >= V.rows() && Fpp(TT[i][2][0], TTi[i][2][0]) >= V.rows()) {//If both vertices are non-original, connect them
					Fp.row(i) << Fpp(i, 2), Fpp(i, 0), Fpp(TT[i][2][0], TTi[i][2][0]);
				}

			}
			else {//Boundary face
				Fp.row(i) = Fpp.row(i);
			}
		}

		Vout = Vp;
		Fout = Fp;

		// Set up the viewer to display the new mesh
		V = Vout; F = Fout;
		viewer.data.clear();
		viewer.data.set_mesh(V, F);
	}

	return true;
}

bool callback_load_mesh(Viewer& viewer, string filename)
{
	igl::readOFF(filename, V, F);
	viewer.data.clear();
	viewer.data.set_mesh(V, F);
	viewer.data.compute_normals();
	viewer.core.align_camera_center(viewer.data.V);
	std::cout << filename.substr(filename.find_last_of("/") + 1) << endl;
	return true;
}

int main(int argc, char *argv[]) {
	// Show the mesh
	Viewer viewer;
	viewer.callback_key_down = callback_key_down;
	//viewer.callback_load_mesh = callback_load_mesh;
	

	if (argc == 2)
	{
		// Read mesh
		igl::readOFF(argv[1], V, F);

	}
	else
	{
		// Read mesh
		igl::readOFF("../data/plane.off", V, F);
	}
	callback_key_down(viewer, '1', 0);

	viewer.data.set_mesh(V, F);
	viewer.core.align_camera_center(V);
	viewer.launch();

}