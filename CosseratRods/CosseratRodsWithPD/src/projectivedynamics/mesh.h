#ifndef __MESH_PD__
#define __MESH_PD__

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "helpers.h"

class Mesh {
public:
	Mesh(const VectorX& variables, const int nDOF) : _variables(variables), _variables_rest(variables), _nDOF(nDOF) {
		_velocity = VectorX::Zero(_nDOF * 3, 1);
		_angular_velocity = VectorX::Zero((_nDOF - 1) * 3, 1);
	}

	int NumDOF() const { return _nDOF; }
	const SparseMatrix& GetMassMatrix() const { return _M; }
	const SparseMatrix& GetMassInverseMatrix() const { return _Minv; }

	void InitTorque(const VectorX& torque) { _torque_rest = torque; SetTorque(torque); }
	void SetTorque(const VectorX& torque) { _torque = torque; }

	virtual void Render() {}

	VectorX GetVariables() const { return _variables; }
	VectorX GetVelocities() const { return _velocity; }
	VectorX GetAngularVelocities() const { return _angular_velocity; }

	void SetVariables(const VectorX& q) { _variables = q; }
	void SetVelocities(const VectorX& v) { _velocity = v; }
	void SetAngularVelocities(const VectorX& w) { _angular_velocity = w; }

protected:
	VectorX _variables, _variables_rest, _velocity, _angular_velocity, _torque, _torque_rest;
	int _nDOF;
	SparseMatrix _M, _Minv;
};

#endif //__MESH_PD__
