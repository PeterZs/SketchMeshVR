#include <utility>
#ifndef __SCENE_PD__
#define __SCENE_PD__

#include "body.h"
#include "helpers.h"
#include "input.h"

class Constraint;

class Scene {

public:
	Body* InitializeScene();

private:
	RodParameters _rod_param;
	SimulationParameters _simulation_param;
	SolverParameters _solver_param;
};

#endif //__SCENE_PD__
