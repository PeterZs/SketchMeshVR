#include <GL/glew.h>
#include "polyline.h"

VectorX RotablePolyLine::generatePolyline(const RodParameters& p) {
	VectorX variables(3 * p.num_samples + 4 * p.num_orientations);
	// Positions
	for (int i = 0; i < p.num_samples; i++) {
		const Scalar x = p.length * ((Scalar)i / (p.num_samples - 1));
		variables.blockVector3(i) = Vector3(x, p.height, (Scalar) 0.0);
	}
	// Quaternions
	for (int i = 0; i < p.num_orientations; i++) {
		const int offset = p.num_samples * 3 + i * 4;
		variables[offset + 0] = (Scalar) 1.0;
		variables[offset + 1] = (Scalar) 0.0;
		variables[offset + 2] = (Scalar) 0.0;
		variables[offset + 3] = (Scalar)  0.0;
		variables.block(offset, 0, 4, 1).normalize();
	}
	return variables;
}

RotablePolyLine::RotablePolyLine(const RodParameters& p, const VectorX& variables, const int& N) : Mesh(variables, N)
{
	SetRadius(p.r);
	buildMasMatrix(variables, p);
}

void RotablePolyLine::Render() {
	glPushAttrib(GL_CURRENT_BIT | GL_LINE_BIT | GL_POINT_BIT);

	/// Render of lines connecting points
	render_lines_polyline();

	/// Render of orientation frames from the quaternions
	render_orientation_frames();

	/// Render the points in the rod
	render_points_polyline();

	glPopAttrib();
}

void RotablePolyLine::buildMasMatrix(const VectorX& variables, const RodParameters& p)
{
	std::vector< Triplet > m, m_inv;
	std::vector<Scalar> weight;

	// mass matrix M for position variables
	for (int k = 0; k < _nDOF; k++) {
		Scalar length_segment = p.length / (_nDOF - 1);
		if (k == 0 || k == _nDOF - 1) {
			length_segment *= 0.5;
		}

		weight.push_back(p.mass_density * M_PI * p.r * p.r * length_segment);
		for (int j = 0; j < 3; j++) {
			m.push_back(Triplet(3 * k + j, 3 * k + j, weight[k]));
			m_inv.push_back(Triplet(3 * k + j, 3 * k + j, 1.0 / weight[k]));
		}
	}

	// inertia matrix J for orientation variables
	for (int k = 0; k < p.num_orientations; k++) {
		Scalar I1, I2, I3, mass;
		mass = (weight[k] + weight[k + 1]) * 0.5;
		I3 = mass * _radius * _radius / (2.0);
		I1 = mass * _radius * _radius / (4.0);
		I2 = I1;

		Scalar I1_inv = 1.0 / I1;
		Scalar I2_inv = 1.0 / I2;
		Scalar I3_inv = 1.0 / I3;

		const int offset = 3 * _nDOF + 4 * k;
		int j = 1;
		m.push_back(Triplet(offset + j, offset + j, I1));
		m_inv.push_back(Triplet(offset + j, offset + j, I1_inv));
		j = 2;
		m.push_back(Triplet(offset + j, offset + j, I2));
		m_inv.push_back(Triplet(offset + j, offset + j, I2_inv));
		j = 3;
		m.push_back(Triplet(offset + j, offset + j, I3));
		m_inv.push_back(Triplet(offset + j, offset + j, I3_inv));
	}

	_M = SparseMatrix(3 * _nDOF + 4 * (_nDOF - 1), 3 * _nDOF + 4 * (_nDOF - 1));
	_Minv = SparseMatrix(3 * _nDOF + 4 * (_nDOF - 1), 3 * _nDOF + 4 * (_nDOF - 1));

	_M.setFromTriplets(m.begin(), m.end());
	_Minv.setFromTriplets(m_inv.begin(), m_inv.end());
}

void RotablePolyLine::render_lines_polyline() {
	glLineWidth(2.0);
	glBegin(GL_LINE_STRIP);
	glColor3f(0.0f, 0.45f, 0.9f);
	for (int i = 0; i < _nDOF; i++) {
		glVertex3dv(_variables.blockVector3(i).data());
	}
	glEnd();
}

void RotablePolyLine::render_orientation_frames() {
	const Scalar lengthFrame = _radius;
	const Vector3 e1(0.0, lengthFrame, 0.0);
	const Vector3 e2(0.0, 0.0, lengthFrame);
	const Vector3 e3(lengthFrame, 0.0, 0.0);
	Quaternion q_e3(0.0, e3[0], e3[1], e3[2]);
	Quaternion q_e1(0.0, e1[0], e1[1], e1[2]);
	Quaternion q_e2(0.0, e2[0], e2[1], e2[2]);

	for (int i = 0; i < _nDOF - 1; i++) {
		const Vector3 origin = (_variables.blockVector3(i) + _variables.blockVector3(i + 1))*0.5;

		Quaternion q(_variables[_nDOF * 3 + i * 4 + 0], _variables[_nDOF * 3 + i * 4 + 1], _variables[_nDOF * 3 + i * 4 + 2], _variables[_nDOF * 3 + i * 4 + 3]);
		const Vector3 d1 = origin + (q * q_e1 * q.inverse()).vec();
		const Vector3 d2 = origin + (q * q_e2 * q.inverse()).vec();
		const Vector3 d3 = origin + (q * q_e3 * q.inverse()).vec();

		glLineWidth(1.4);
		glBegin(GL_LINE_STRIP);
		glColor3f(0.0f, 1.0f, 1.0f);
		glVertex3dv(origin.data());
		glVertex3dv(d3.data());
		glEnd();

		glLineWidth(1.4);
		glBegin(GL_LINE_STRIP);
		glColor3f(1.0f, 0.0f, 0.0f);
		glVertex3dv(origin.data());
		glVertex3dv(d1.data());
		glEnd();

		glLineWidth(1.4);
		glBegin(GL_LINE_STRIP);
		glColor3f(0.f, 0.7f, 0.0f);
		glVertex3dv(origin.data());
		glVertex3dv(d2.data());
		glEnd();
	}
}

void RotablePolyLine::render_points_polyline() {
	glPointSize(3.25);
	glBegin(GL_POINTS);
	glColor3f(0.f, 0.f, 0.5f);
	for (int i = 0; i < _nDOF; i++)
		glVertex3dv(_variables.blockVector3(i).data());
	glEnd();
}

