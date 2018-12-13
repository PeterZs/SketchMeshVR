#include "softBody.h"
#include "input.h"

void SoftBody::SetSimulationParameters(const SolverParameters& solver_p, const SimulationParameters& sim_param) {
	_solver_param = solver_p;
	_simulation_param = sim_param;
	_h = _simulation_param.h;

	/// init linear external forces
	_f_ext = VectorX(_mesh->NumDOF() * 3);
	for (int i = 0; i < _mesh->NumDOF(); i++) {
		_f_ext.blockVector3(i) = _mesh->GetMassMatrix().block(i, i, 3, 3) * _simulation_param.gravity;
	}

	/// init rotational external forces
	_torque = VectorX(_mesh->NumDOF() * 3);
	for (int i = 0; i < _mesh->NumDOF(); i++) {
		_torque.blockVector3(i) = _simulation_param.torque;
	}
	_mesh->InitTorque(_torque);
}

void SoftBody::Render() {
	_mesh->Render();
}

void SoftBody::InitLHS() {
	computeLHS(_constraints);
}

void SoftBody::computeLHS(const std::vector<Constraint*>& constraints) {
	SparseMatrix M = _mesh->GetMassMatrix() / (_h * _h);
	for (int i = 0; i < constraints.size(); i++) {
		const SparseMatrix A = constraints[i]->MatrixA(); //m x m
		const SparseMatrix S = constraints[i]->MatrixS(); //m x n
		const Scalar w = constraints[i]->Stiffness();
		M += w * (S.transpose() * A.transpose() * A * S);
	}
	_LHS.compute(M);
}

void SoftBody::AddToRightHandSide(const std::vector<Constraint*>& constraints, const VectorX q, VectorX& RHS) const {
	for (int j = 0; j < constraints.size(); j++) { //Equation (10), right
		const VectorX p = constraints[j]->Evaluate(q);
		const SparseMatrix& SAB = constraints[j]->MatrixSAB();
		RHS += SAB * p;
	}
}

void SoftBody::AddConstraint(Constraint* c) {
	c->computeSAB();
	_constraints.push_back(c);
}


std::pair<VectorX, VectorX> obtainVariablesFromVector(const VectorX& variables) {
	int N = (variables.size() + 4) / 7;

	VectorX x0(N * 3);
	for (int i = 0; i < N; i++) {
		x0.blockVector3(i) = variables.blockVector3(i);
	}
	VectorX u0 = VectorX((N - 1) * 4);
	for (int i = 0; i < N - 1; i++) {
		const int offset = N * 3 + i * 4;
		u0.blockVector4(i) = Vector4(variables[offset], variables[offset + 1], variables[offset + 2], variables[offset + 3]);
		u0.blockVector4(i).normalize();
	}
	return { x0, u0 };
}

VectorX buildVectorWithVariables(const VectorX& x, const VectorX& u) {
	VectorX s0(x.size() + u.size());
	for (int i = 0; i < x.size() / 3; i++) {
		s0.blockVector3(i) = x.blockVector3(i);
	}
	for (int i = 0; i < u.size() / 4; i++) {
		s0.block(x.size() + i * 4, 0, 4, 1) = u.blockVector4(i);
	}
	return s0;
}

void SoftBody::Integrate() {
	const auto M = _mesh->GetMassMatrix();

	/// returns a pair with {position, orientation} variables
	auto variables = obtainVariablesFromVector(_mesh->GetVariables());

	/// position variables update
	auto x1 = updateLinearMomentum(variables.first);
	/// orientation variables update
	auto u1 = updateAngularMomentum(variables.second);
	/// Initialize s0: system variables with x1 and u1
	auto s0 = buildVectorWithVariables(x1, u1);

	VectorX q1 = s0; /// Initialize q1 with explicit integration step. 
	for (int i = 0; i < _solver_param.num_iterations; ++i) { ///q1 is iteratively updated with the local projections.
		VectorX RHS = (M / (_h * _h)) * s0;
		AddToRightHandSide(_constraints, q1, RHS);
		q1 = _LHS.solve(RHS);
	}

	/// Separate system variables q1 into positions and orientations vectors
	auto solution_solver = obtainVariablesFromVector(q1);
	auto updated_velocities = updateVelocities(solution_solver, variables);

	_mesh->SetVariables(q1);								/// Update mesh variables with the converged system variables q1
	_mesh->SetVelocities(updated_velocities.first);			/// update linear velocity
	_mesh->SetAngularVelocities(updated_velocities.second); /// update angular velocity
}


VectorX SoftBody::updateLinearMomentum(const VectorX& x0) {
	const int N = _mesh->NumDOF();

	const auto mass_matrix_inv = _mesh->GetMassInverseMatrix();
	const auto M_inv = mass_matrix_inv.block(0, 0, N * 3, N * 3);

	const auto v0 = _mesh->GetVelocities();			  // v_{t}
	return x0 + _h * v0 + _h * _h * M_inv * _f_ext;   // v_{t+1}
}

VectorX SoftBody::updateAngularMomentum(const VectorX& u0) {
	const int N = _mesh->NumDOF();
	const VectorX w0 = _mesh->GetAngularVelocities(); //w_{t}
	const auto J = _mesh->GetMassMatrix().block(N * 3, N * 3, (N - 1) * 4, (N - 1) * 4);
	const auto J_inv = _mesh->GetMassInverseMatrix().block(N * 3, N * 3, (N - 1) * 4, (N - 1) * 4);
	VectorX u1((N - 1) * 4);

	for (int i = 0; i < N - 1; i++) {

		auto u_t = u0.blockVector4(i);			// orientation
		auto torque = _torque.blockVector3(i);	// torque
		auto w_t = w0.blockVector3(i);			// angular velocity

		const Matrix3 Jn = J.block(i * 4 + 1, i * 4 + 1, 3, 3);
		const Matrix3 Jn_inv = J_inv.block(i * 4 + 1, i * 4 + 1, 3, 3);

		const Vector3 s_w_t = w_t + _h * Jn_inv * (torque - w_t.cross(Jn * w_t));
		const Quaternion s_w_t_q(0.0, s_w_t[0], s_w_t[1], s_w_t[2]);

		const Quaternion u_t_q(u_t[0], u_t[1], u_t[2], u_t[3]);
		const Quaternion product = u_t_q * s_w_t_q;
		const Vector4 product_vec(product.w(), product.vec()[0], product.vec()[1], product.vec()[2]);

		u1.blockVector4(i) = u_t + 0.5 * _h * product_vec;
		u1.blockVector4(i).normalize();
	}
	return u1;
}

std::pair<VectorX, VectorX> SoftBody::updateVelocities(const std::pair<VectorX, VectorX>& solution, const std::pair<VectorX, VectorX>& previous_value) {
	auto x1 = solution.first;
	auto u1 = solution.second;
	auto x0 = previous_value.first;
	auto u0 = previous_value.second;

	VectorX v1 = (x1 - x0) / _h;

	VectorX w1(3 * (_mesh->NumDOF() - 1));
	for (int i = 0; i < _mesh->NumDOF() - 1; i++) {
		Quaternion u_0(u0.blockVector4(i)[0], u0.blockVector4(i)[1], u0.blockVector4(i)[2], u0.blockVector4(i)[3]);
		Quaternion u_1(u1.blockVector4(i)[0], u1.blockVector4(i)[1], u1.blockVector4(i)[2], u1.blockVector4(i)[3]);
		Quaternion product = u_0.inverse() * u_1;
		w1.blockVector3(i) = product.vec() * 2.0 / _h;
	}
	return { v1, w1 };
}