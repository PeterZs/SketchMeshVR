#include "helpers.h"
#include "constraint.h"

AttachmentConstraint::AttachmentConstraint(const int& idx, const VectorX& q, const int& stiffness) : Constraint(stiffness) {
	_fixed_variable = q.blockVector3(idx);
	_idx = 3 * idx;

	std::vector< Triplet > atriplet = {
		{ 0, 0, 1.0 },
		{ 1, 1, 1.0 },
		{ 2, 2, 1.0 }
	};
	_Ai = SparseMatrix(3, 3);
	_Ai.setFromTriplets(atriplet.begin(), atriplet.end());

	std::vector< Triplet > btriplet = {
		{ 0, 0, 1.0 },
		{ 1, 1, 1.0 },
		{ 2, 2, 1.0 },
	};
	_Bi = SparseMatrix(3, 3);
	_Bi.setFromTriplets(btriplet.begin(), btriplet.end());

	std::vector< Triplet > striplet = {
		{ 0, idx * 3 + 0, 1.0 },
		{ 1, idx * 3 + 1, 1.0 },
		{ 2, idx * 3 + 2, 1.0 }
	};
	_Si = SparseMatrix(3, q.size());
	_Si.setFromTriplets(striplet.begin(), striplet.end());
}

VectorX AttachmentConstraint::Evaluate(const VectorX& q) {
	return _fixed_variable;
}

StretchShearConstraint::StretchShearConstraint(const int& a, const int& b, const VectorX& q, const RodParameters& rod_parameters) : Constraint(ComputeStretchStiffness(a, b, q, rod_parameters)), _a(a), _b(b) {
	auto length_rest = (q.blockVector3(a) - q.blockVector3(b)).norm();
	_nDOF = rod_parameters.num_samples;
	_orientationID = _nDOF * 3 + a * 4;

	const Scalar cte_a = -1.0 / length_rest;
	const Scalar cte_b = 1.0 / length_rest;
	const Scalar cte_u = 1.0;
	std::vector< Triplet > atriplet = {
		{ 0, 0, cte_a },
		{ 1, 1, cte_a },
		{ 2, 2, cte_a },
		{ 0, 3, cte_b },
		{ 1, 4, cte_b },
		{ 2, 5, cte_b },
		{ 3, 6, cte_u },
		{ 4, 7, cte_u },
		{ 5, 8, cte_u },
		{ 6, 9, cte_u }
	};
	_Ai = SparseMatrix(7, 10);
	_Ai.setFromTriplets(atriplet.begin(), atriplet.end());

	std::vector< Triplet > btriplet = {
		{ 0, 0, 1.0 },
		{ 1, 1, 1.0 },
		{ 2, 2, 1.0 },
		{ 3, 3, cte_u },
		{ 4, 4, cte_u },
		{ 5, 5, cte_u },
		{ 6, 6, cte_u }
	};
	_Bi = SparseMatrix(7, 7);
	_Bi.setFromTriplets(btriplet.begin(), btriplet.end());
	std::vector< Triplet > striplet = {
		{ 0, a * 3 + 0, 1.0 },
		{ 1, a * 3 + 1, 1.0 },
		{ 2, a * 3 + 2, 1.0 },
		{ 3, b * 3 + 0, 1.0 },
		{ 4, b * 3 + 1, 1.0 },
		{ 5, b * 3 + 2, 1.0 },
		{ 6, _orientationID + 0, cte_u },
		{ 7, _orientationID + 1, cte_u },
		{ 8, _orientationID + 2, cte_u },
		{ 9, _orientationID + 3, cte_u }
	};
	_Si = SparseMatrix(10, q.size());
	_Si.setFromTriplets(striplet.begin(), striplet.end());
}

Scalar StretchShearConstraint::ComputeStretchStiffness(const int& a, const int& b, const VectorX& q, const RodParameters& rod_parameters) {
	Scalar length_rest = (q.blockVector3(a) - q.blockVector3(b)).norm();
	Scalar E = rod_parameters.E;
	Scalar r = rod_parameters.r;

	Scalar stiffness = E * M_PI * r * r * length_rest;
	return stiffness;
}

VectorX StretchShearConstraint::Evaluate(const VectorX& q) {

	Vector3 diff = (q.blockVector3(_b) - q.blockVector3(_a)).normalized();
	const Vector4 quat_vec(q[_orientationID], q[_orientationID + 1], q[_orientationID + 2], q[_orientationID + 3]);
	Quaternion e3(0.0, 1.0, 0.0, 0.0);
	Quaternion quat(q[_orientationID], q[_orientationID + 1], q[_orientationID + 2], q[_orientationID + 3]);
	Vector3 d3 = (quat * e3 * quat.inverse()).vec();

	VectorX p(7);
	p.blockVector3(0) = d3;

	const  Quaternion diff_q = Quaternion::FromTwoVectors(d3.normalized(), diff);
	const Quaternion quat_updated = diff_q * quat;
	const Vector4 vec_quat_updated(quat_updated.w(), quat_updated.vec()[0], quat_updated.vec()[1], quat_updated.vec()[2]);
	p.block(3, 0, 4, 1) = vec_quat_updated.normalized();

	return p;
}

BendTwistConstraint::BendTwistConstraint(const int& a, const VectorX& q, const RodParameters& rod_parameters) : Constraint(ComputeBendStiffness(a, q, rod_parameters)), _a(a) {

	_nDOF = rod_parameters.num_samples;
	_q_idx = _nDOF * 3 + (_a - 1) * 4;
	_u_idx = _nDOF * 3 + _a * 4;

	_q0 = Quaternion(q[_q_idx], q[_q_idx + 1], q[_q_idx + 2], q[_q_idx + 3]);
	_u0 = Quaternion(q[_u_idx], q[_u_idx + 1], q[_u_idx + 2], q[_u_idx + 3]);

	std::vector< Triplet > atriplet = {
		{ 0, 0, 1.0 },
		{ 1, 1, 1.0 },
		{ 2, 2, 1.0 },
		{ 3, 3, 1.0 },
		{ 4, 4, 1.0 },
		{ 5, 5, 1.0 },
		{ 6, 6, 1.0 },
		{ 7, 7, 1.0 }
	};
	_Ai = SparseMatrix(8, 8);
	_Ai.setFromTriplets(atriplet.begin(), atriplet.end());

	std::vector< Triplet > btriplet = {
		{ 0, 0, 1.0 },
		{ 1, 1, 1.0 },
		{ 2, 2, 1.0 },
		{ 3, 3, 1.0 },
		{ 4, 4, 1.0 },
		{ 5, 5, 1.0 },
		{ 6, 6, 1.0 },
		{ 7, 7, 1.0 }
	};
	_Bi = SparseMatrix(8, 8);
	_Bi.setFromTriplets(btriplet.begin(), btriplet.end());

	std::vector< Triplet > striplet = {
		{ 0, _q_idx + 0, 1.0 },
		{ 1, _q_idx + 1, 1.0 },
		{ 2, _q_idx + 2, 1.0 },
		{ 3, _q_idx + 3, 1.0 },
		{ 4, _u_idx + 0, 1.0 },
		{ 5, _u_idx + 1, 1.0 },
		{ 6, _u_idx + 2, 1.0 },
		{ 7, _u_idx + 3, 1.0 }
	};
	_Si = SparseMatrix(8, q.size());
	_Si.setFromTriplets(striplet.begin(), striplet.end());
}

Scalar BendTwistConstraint::ComputeBendStiffness(const int& a, const VectorX& q, const RodParameters& rod_parameters) {
	Scalar length_rest = (q.blockVector3(a - 1) - q.blockVector3(a)).norm();
	Scalar E = rod_parameters.E;
	Scalar vau = rod_parameters.vau;
	Scalar r = rod_parameters.r;

	Scalar stiffness = (E / (1.0 + vau)) * M_PI * r * r * r * r / length_rest;
	return stiffness;
}

VectorX BendTwistConstraint::Evaluate(const VectorX& q) {

	Quaternion quat_q(q[_q_idx], q[_q_idx + 1], q[_q_idx + 2], q[_q_idx + 3]);
	const Quaternion quat_u(q[_u_idx], q[_u_idx + 1], q[_u_idx + 2], q[_u_idx + 3]);

	const Quaternion omega_quat0 = _q0.conjugate() * _u0;

	const Vector3 omega_vec = (omega_quat0.inverse() * quat_q.conjugate() * quat_u).vec() * 0.5;
	const Vector4 omega_vec4 = (Vector4(1.0, omega_vec[0], omega_vec[1], omega_vec[2])).normalized();

	const Quaternion omega_quat(omega_vec4[0], omega_vec4[1], omega_vec4[2], omega_vec4[3]);
	const Quaternion quat_q_update = quat_q * omega_quat;
	const Quaternion quat_u_update = quat_u * omega_quat.inverse();

	Vector4 vec_q_update(quat_q_update.w(), quat_q_update.vec()[0], quat_q_update.vec()[1], quat_q_update.vec()[2]);
	Vector4 vec_u_update(quat_u_update.w(), quat_u_update.vec()[0], quat_u_update.vec()[1], quat_u_update.vec()[2]);

	VectorX p(8);
	p.blockVector4(0) = vec_q_update.normalized();
	p.blockVector4(1) = vec_u_update.normalized();
	return p;
}
