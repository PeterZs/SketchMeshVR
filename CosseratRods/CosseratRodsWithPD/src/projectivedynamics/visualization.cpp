#include "visualization.h"
#include "helpersrendering.h"

int screenWidth = 1280, screenHeight = 720, cur_button = -1, last_x, last_y;;
bool pressedCTRL = false, pressedSHIFT = false;
Vector3 eye, lookat, up;
Matrix4 P;
bool simulationPaused = true;
Body* _body;
int timesteps_passed = 0;

void Anim() {
	if (!simulationPaused) {
		_body->Integrate();
		timesteps_passed++;
	}
	glutPostRedisplay();
}

void initGL() {
	glClearColor(0.95f, 0.95f, 0.95f, 0.0);
	glEnable(GL_DEPTH_TEST);

	glEnable(GL_LIGHT0);
	glEnable(GL_COLOR_MATERIAL);
}

void Display() {
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glMultMatrixd(P.data());
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(eye[0], eye[1], eye[2], lookat[0], lookat[1], lookat[2], up[0], up[1], up[2]);

	RenderGrid<Scalar>(1.0f);

	_body->Render();
	glutSwapBuffers();
}

void Reshape(int w, int h) {
	screenWidth = w;
	screenHeight = h;
	glViewport(0, 0, screenWidth, screenHeight);
	P = Perspective<Scalar>(30.0, (GLdouble)screenWidth / screenHeight, 0.1, 20.0);
}

void MouseMotion(int x, int y) {
	int dx = x - last_x;
	int dy = y - last_y;
	switch (cur_button) {
	case GLUT_LEFT_BUTTON:
		if (pressedCTRL) {
			eye = RotateCamera<Scalar>(eye, lookat, dx, dy);
		}
		else if (pressedSHIFT) {
			std::tie(eye, lookat) = MoveCamera<Scalar>(eye, lookat, up, P, Eigen::Vector2i(x, y), Eigen::Vector2i(last_x, last_y), screenWidth, screenHeight);
		}
		else {
			eye = ZoomCamera<Scalar>(eye, lookat, dx);
		}
		break;
	}
	last_x = x;
	last_y = y;
	glutPostRedisplay();
}


void Keyboard(unsigned char key, int x, int y) {
	Scalar displacement = 1.0e-3;
	switch (key) {
	case 27:
		exit(0);
	case ' ':
		simulationPaused = !simulationPaused;
		PRINTVAR(timesteps_passed);
		break;
	}
	glutPostRedisplay();
}

void Mouse(int button, int state, int x, int y) {
	cur_button = button;
	pressedCTRL = glutGetModifiers() & GLUT_ACTIVE_CTRL ? true : false;
	pressedSHIFT = glutGetModifiers() & GLUT_ACTIVE_SHIFT ? true : false;
	if (state == GLUT_DOWN) {
		last_x = x;
		last_y = y;
	}
}


int initialize_renderer(int argc, char **argv, Body* body) {
	_body = body;
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
	glutInitWindowSize(screenWidth, screenHeight);
	glutInitWindowPosition(0, 0);
	glutCreateWindow("Cosserat Rods with PD");

	initGL();

	eye = Vector3(2.0, 2.0, 2.0);
	lookat = Vector3(0.0, 0.0, 0.0);
	up = Vector3(0.0, 1.0, 0.0);

	glutDisplayFunc(Display);
	glutMotionFunc(MouseMotion);
	glutMouseFunc(Mouse);
	glutReshapeFunc(Reshape);
	glutKeyboardFunc(Keyboard);
	glutIdleFunc(Anim);
	glutMainLoop();

	return 0;
}
