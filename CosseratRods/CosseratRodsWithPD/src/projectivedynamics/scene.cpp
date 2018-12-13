#include "scene.h"
#include "polyline.h"
#include "softBody.h"
#include "mesh.h"
#include "constraint.h"

Body* Scene::InitializeScene()
{
	auto variables = RotablePolyLine::generatePolyline(_rod_param);
	int N = (variables.size() + 4) / 7;
	RotablePolyLine* polyline = new RotablePolyLine(_rod_param, variables, N);

	SoftBody* softBody = new SoftBody(polyline);
	softBody->SetSimulationParameters(_solver_param, _simulation_param);

	/// Attachment constraints
	if (_rod_param.attachment.first) softBody->AddConstraint(new AttachmentConstraint(0, variables, 1.0e5));
	if (_rod_param.attachment.second) softBody->AddConstraint(new AttachmentConstraint(_rod_param.num_samples - 1, variables, 1.0e5));

	/// Cosserat Stretch and Shear constraints
	for (int i = 0; i < _rod_param.num_samples - 1; i++) {
		softBody->AddConstraint(new StretchShearConstraint(i, i + 1, variables, _rod_param));
	}

	/// Cosserat Bend and Twist constraints
	for (int i = 1; i < _rod_param.num_orientations; i++) {
		softBody->AddConstraint(new BendTwistConstraint(i, variables, _rod_param));
	}

	softBody->InitLHS();

	return softBody;
}

