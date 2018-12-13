#include <utility>
#include "helpers.h"

#ifndef __INPUT_PD__
#define __INPUT_PD__

struct RodParameters {
	Scalar E{ (Scalar) 1.0e5 };				/// Young's Modulus (kPa = g m-1 s-2) [1e5 kPa = rubber, 7e7 kPa = aluminium]
	Scalar r{ (Scalar) 0.0010 };				/// Rod's radius (m)
	Scalar mass_density{ (Scalar)1.5e3 };	/// mass density (kg m-3) Given in kg because forces handle kgs and not grams [1.5e3 = rubber, 2.7e3 = aluminium (decrease h to 0.005 for alu)]
	Scalar vau{ (Scalar) 0.48 };			/// Poisson ratio [0 = cork, 0.3 = metal, 0.48 = rubber]

	int num_samples{ 20 };
	int num_orientations{ num_samples - 1 };
	Scalar length{ (Scalar) 1.0 };			/// length of the rod (m)
	Scalar height{ (Scalar) 1.0 };			/// height of initialized rod (m)
	std::pair<bool, bool> attachment{ true, true }; /// indicates whether the endpoint is attached {left endpoint, right endpoint}
};
struct SimulationParameters {
	Vector3 gravity{ (Scalar) 0.0,  (Scalar)-9.81, (Scalar) 0.0 };		/// linear external force (m s-2)
	Vector3 torque = Vector3::Zero();		/// rotational external force
	Scalar h{ (Scalar) 0.01 };				/// timestep (s)
	//Scalar friction{ (Scalar) 8.4 };		/// friction
};
struct SolverParameters {
	int num_iterations{ 20 };		/// solver iterations
};
#endif //__INPUT_PD__
