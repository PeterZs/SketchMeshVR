#ifndef __CONSTRAINT_PD_
#define __CONSTRAINT_PD_

#include <Eigen/Sparse>
#include "input.h"

class Constraint {
	/*Class containing the different types of constraints
		- See [1]: Projective Dynamics: Fusing Constraint Projections for Fast Simulation, Bouaziz et al. (https://www.cs.utah.edu/~ladislav/bouaziz14projective/bouaziz14projective.pdf) */
public:
	Constraint(const Scalar stiffness) : _stiffness(stiffness) {}

	// Return the projection variables p_i in the local solve -- see [1], eq.9
	virtual VectorX Evaluate(const VectorX& q) { return  q; }

	// A,B,S are matrices for setting up the global solve -- see [1], eq.10
	// A,B,S are defined per constraint (each constraint is i in eq.10)
	const SparseMatrix& MatrixA() const { return _Ai; }
	const SparseMatrix& MatrixB() const { return _Bi; }
	const SparseMatrix& MatrixS() const { return _Si; }
	const SparseMatrix& MatrixSAB() const { return _SABi; }

	// Stiffness is the weight w_i of the potential in the global solve -- see[1], eq. 10
	virtual Scalar Stiffness() const = 0;

	//	Precompute the matrix product to reduce computation time
	void computeSAB() { _SABi = _stiffness * _Si.transpose() * _Ai.transpose() * _Bi; }

protected:
	SparseMatrix _Ai, _Bi, _Si, _SABi; //A, B as defined in Bouaziz14
	Scalar _stiffness;
};

/// Attachment Constraint
class AttachmentConstraint : public Constraint {
public:
	AttachmentConstraint(const int& idx, const VectorX& q, const int& stiffness);
	VectorX Evaluate(const VectorX& q) override;
	Scalar Stiffness() const override { return _stiffness; }

private:
	int _idx;
	VectorX _fixed_variable;
};

/// Cosserat stretch and shear constraint
class StretchShearConstraint : public Constraint {
public:
	StretchShearConstraint(const int& a, const int& b, const VectorX& q, const RodParameters& rod_parameters);
	VectorX Evaluate(const VectorX& q) override;
	Scalar ComputeStretchStiffness(const int& a, const int& b, const VectorX& q, const RodParameters& rod_parameters);
	Scalar Stiffness() const override { return _stiffness; }

private:
	int _a, _b, _orientationID, _nDOF;
};

/// Cosserat stretch and shear constraint
class BendTwistConstraint : public Constraint {
public:
	BendTwistConstraint(const int& a, const VectorX& q, const RodParameters& rod_parameters);
	VectorX Evaluate(const VectorX& q) override;
	Scalar ComputeBendStiffness(const int& a, const VectorX& q, const RodParameters& rod_parameters);
	Scalar Stiffness() const override { return _stiffness; }

private:
	int _a, _q_idx, _u_idx, _nDOF;
	Quaternion _q0, _u0;
};

#endif //__CONSTRAINT_PD_
