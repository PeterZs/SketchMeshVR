#include "helpers.h"
#include <GL/glew.h>
#include <freeglut.h>

template<typename T, size_t D>
struct Ray {
	Ray(const Eigen::Matrix<T, D, 1>& origin, const Eigen::Matrix<T, D, 1>& direction) : origin(origin), direction(direction) {}
	Eigen::Matrix<T, D, 1> origin, direction;
};

template<typename T, size_t D>
struct Plane {
	Plane(const Eigen::Matrix<T, D, 1>& p, const Eigen::Matrix<T, D, 1>& n) : p(p), n(n) {}
	Eigen::Matrix<T, D, 1> p, n;
};

template<class T>
Eigen::Matrix<T, 4, 4> Frustum(const T& l, const T& r, const T& b, const T& t, const T& n, const T& f) {
	Eigen::Matrix<T, 4, 4> F;
	F << (T)2 * n / (r - l), (T)0, (r + l) / (r - l), (T)0,
		(T)0, ((T)2 * n) / (t - b), (t + b) / (t - b), (T)0,
		(T)0, (T)0, -(f + n) / (f - n), (-(T)2 * f * n) / (f - n),
		(T)0, (T)0, (T)-1, (T)0;
	return F;
}

template<class T>
Eigen::Matrix<T, 4, 4> Perspective(T fovy, T aspect, T n, T f) {
	const T pi = (T)3.14159265359;
	T xmin, xmax, ymin, ymax;
	ymax = n * std::tan(fovy * pi / (T)360);
	ymin = -ymax;
	xmin = ymin * aspect;
	xmax = ymax * aspect;
	return Frustum<T>(xmin, xmax, ymin, ymax, n, f);
}

template<class T>
void RenderGrid(const T scale = (T)10, const int n = 21) {
	glPushAttrib(GL_ENABLE_BIT | GL_CURRENT_BIT | GL_LINE_BIT);
	glDisable(GL_LIGHTING);
	glDisable(GL_TEXTURE_2D);
	// 			glDisable(GL_DEPTH_TEST);
	glLineWidth(1.0);
	glColor3f(192.0f / 255.0f, 192.0f / 255.0f, 192.0f / 255.0f);
	glBegin(GL_LINES);
	for (int i = 0; i < n; i++) {
		const T t = (T)i / (n - 1);
		const T f = (T)2 * scale * t - scale;
		if (i != (n - 1) / 2) {
			glVertex3f((float)f, 0.0f, (float)-scale);
			glVertex3f((float)f, 0.0f, (float)scale);
			glVertex3f((float)-scale, 0.0f, (float)f);
			glVertex3f((float)scale, 0.0f, (float)f);
		}
	}
	glEnd();
	glColor3f(140.0f / 255.0f, 140.0f / 255.0f, 140.0f / 255.0f);
	glBegin(GL_LINES);
	glVertex3f(0.0f, 0.0f, (float)-scale);
	glVertex3f(0.0f, 0.0f, (float)scale);
	glVertex3f((float)-scale, 0.0f, 0.0f);
	glVertex3f((float)scale, 0.0f, 0.0f);
	glEnd();
	glPopAttrib();
}

template<class T>
Eigen::Matrix<T, 4, 4> LookAt(const Eigen::Matrix<T, 3, 1>& pos, const Eigen::Matrix<T, 3, 1>& center, const Eigen::Matrix<T, 3, 1>& up) {
	Eigen::Matrix<T, 3, 1> f = center - pos;
	f.normalize();

	Eigen::Matrix<T, 3, 1> s = f.cross(up);
	s.normalize();

	Eigen::Matrix<T, 3, 1> u = s.cross(f);

	Eigen::Matrix<T, 4, 4> C;
	C << s[0], s[1], s[2], -s[0] * pos.x() - s[1] * pos.y() - s[2] * pos.z(),
		u[0], u[1], u[2], -u[0] * pos.x() - u[1] * pos.y() - u[2] * pos.z(),
		-f[0], -f[1], -f[2], f[0] * pos.x() + f[1] * pos.y() + f[2] * pos.z(),
		0.0f, 0.0f, 0.0f, 1.0f;

	return C;
}

template<class T>
Eigen::Matrix<T, 3, 1> RotateCamera(const Eigen::Matrix<T, 3, 1>& eye, const Eigen::Matrix<T, 3, 1>& lookat, const int dx, const int dy) {
	const Eigen::Matrix<T, 3, 1> neye = eye - lookat;
	const T thetax = -dx * (T)0.007;
	Eigen::Matrix<T, 3, 1> f = -Eigen::Matrix<T, 3, 1>(std::cos(thetax) * neye[0] + std::sin(thetax) * neye[2], neye[1], -std::sin(thetax) * neye[0] + std::cos(thetax) * neye[2]);
	Eigen::Matrix<T, 3, 1> u = Eigen::Matrix<T, 3, 1>((T)0, (T)1, (T)0);
	Eigen::Matrix<T, 3, 1> r = f.cross(u);
	u = r.cross(f);
	const T len = f.norm();
	f.normalize();
	u.normalize();
	const T thetay = -dy * (T)0.007;
	return lookat - len * (std::cos(thetay) * f + std::sin(thetay) * u);
}

template<class T>
std::pair<Eigen::Matrix<T, 3, 1>, Eigen::Matrix<T, 3, 1> > MakeRay(const Eigen::Matrix<T, 2, 1>& pixel, const Eigen::Matrix<T, 4, 4>& Pinv, const Eigen::Matrix<T, 4, 4>& Cinv, const int width, const int height) {
	//pixel in relative image coordinates
	const Eigen::Matrix<T, 2, 1> rel = Eigen::Matrix<T, 2, 1>((T)pixel[0] / width + (T)1 / ((T)2 * width), (T)(height - 1 - pixel[1]) / height + (T)1 / ((T)2 * height));

	//pixel in normalized device coordinates
	const Eigen::Matrix<T, 4, 1> prj = Eigen::Matrix<T, 4, 1>((T)2 * rel[0] - (T)1, (T)2 * rel[1] - (T)1, -(T)1, (T)1);

	//pixel in eye space
	const Eigen::Matrix<T, 4, 1> eye = Pinv * prj;

	//pixel in world space
	const Eigen::Matrix<T, 4, 1> world = Cinv * (eye / eye[3]);

	//ray origin
	const Eigen::Matrix<T, 4, 1> origin = Cinv * Eigen::Matrix<T, 4, 1>((T)0, (T)0, (T)0, (T)1);

	//ray direction
	Eigen::Matrix<T, 3, 1> direction = Eigen::Matrix<T, 3, 1>(world[0], world[1], world[2]) - Eigen::Matrix<T, 3, 1>(origin[0], origin[1], origin[2]);
	direction.normalize();

	return std::make_pair(Eigen::Matrix<T, 3, 1>(origin[0], origin[1], origin[2]), direction);
}

template<class T, size_t D>
bool Intersect(const Ray<T, D>& r, const Plane<T, D>& p, T& t) {
	const T distance = p.n.dot(p.p);
	const T denom = p.n.dot(r.direction);
	if (std::abs(denom) > (T)1.e-6) {
		const T t0 = (distance - p.n.dot(r.origin)) / denom;
		if (t0 >= 0.0f) {
			t = t0;
			return true;
		}
	}
	return false;
}

template<class T>
std::pair<Eigen::Matrix<T, 3, 1>, Eigen::Matrix<T, 3, 1> > MoveCamera(const Eigen::Matrix<T, 3, 1>& eye, const Eigen::Matrix<T, 3, 1>& lookat, const Eigen::Matrix<T, 3, 1>& up, const Eigen::Matrix<T, 4, 4>& P, const Eigen::Vector2i& screenPosition, const Eigen::Vector2i& previousScreenPosition, const int screenWidth, const int screenHeight) {
	const Eigen::Matrix<T, 4, 4> Pinv = P.inverse();
	const Eigen::Matrix<T, 4, 4> Cinv = LookAt<T>(eye, lookat, up).inverse();
	const std::pair< Eigen::Matrix<T, 3, 1>, Eigen::Matrix<T, 3, 1> > ray0 = MakeRay<T>(screenPosition.cast<T>(), Pinv, Cinv, screenWidth, screenHeight);
	const std::pair< Eigen::Matrix<T, 3, 1>, Eigen::Matrix<T, 3, 1> > ray1 = MakeRay<T>(previousScreenPosition.cast<T>(), Pinv, Cinv, screenWidth, screenHeight);

	T t0, t1;
	Eigen::Matrix<T, 3, 1> d((T)0, (T)0, (T)0);
	const Plane<T, 3> xyPlane(Eigen::Matrix<T, 3, 1>((T)0, (T)0, (T)0), Eigen::Matrix<T, 3, 1>((T)0, (T)1, (T)0));
	if (Intersect(Ray<T, 3>(ray0.first, ray0.second), xyPlane, t0) && Intersect(Ray<T, 3>(ray1.first, ray1.second), xyPlane, t1)) {
		const Eigen::Matrix<T, 3, 1> p0 = ray0.first + t0 * ray0.second;
		const Eigen::Matrix<T, 3, 1> p1 = ray1.first + t1 * ray1.second;
		d = p1 - p0;
	}

	return std::make_pair(eye + d, lookat + d);
}

template<class T>
Eigen::Matrix<T, 3, 1> ZoomCamera(const Eigen::Matrix<T, 3, 1>& eye, const Eigen::Matrix<T, 3, 1>& lookat, const int dx) {
	Eigen::Matrix<T, 3, 1> f = lookat - eye;
	T len = f.norm();
	len -= std::sqrt(len) * dx * (T)0.03;
	if (len < 0.1)
		return eye;
	f.normalize();
	return lookat - len * f;
}
