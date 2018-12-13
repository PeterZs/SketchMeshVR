#ifndef __SOFTBODY_PD__
#define __SOFTBODY_PD__

#include "body.h"
#include "mesh.h"
#include "scene.h"
#include "constraint.h"
#include <vector>

class SoftBody : public Body {
public:
	SoftBody(Mesh* mesh) : _mesh(mesh) {}
	~SoftBody() {}
	void Render() override;
	void SetSimulationParameters(const SolverParameters& solver_p, const SimulationParameters& sim_param);

	void Integrate() override;
	void InitLHS();
	void AddConstraint(Constraint* c);

private:
	void computeLHS(const std::vector<Constraint*>& constraints);
	void AddToRightHandSide(const std::vector<Constraint*>& constraints, const VectorX q, VectorX& RHS) const;
	VectorX updateLinearMomentum(const VectorX& x0);
	VectorX updateAngularMomentum(const VectorX& u0);
	std::pair<VectorX, VectorX> updateVelocities(const std::pair<VectorX, VectorX>& solution, const std::pair<VectorX, VectorX>& previous_value);

	std::vector<Constraint*> _constraints;
	LLT _LHS;

	SimulationParameters _simulation_param;
	SolverParameters _solver_param;
	Mesh* _mesh;
	VectorX _f_ext;
	VectorX _torque;
	Scalar _h;

	int timesteps_passed = 0;
};
#endif //__SOFTBODY_PD__