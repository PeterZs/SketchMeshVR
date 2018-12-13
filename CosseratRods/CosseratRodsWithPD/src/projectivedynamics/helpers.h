#ifndef __HELPERS_PD__
#define __HELPERS_PD__

#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Sparse>

#define _USE_MATH_DEFINES
#include <math.h>

#define PRINTVAR(x) std::cout << #x << ": " << x << std::endl
#define PRINTVEC(x) std::cout << #x << ": " << x.transpose() << std::endl
#define PRINTMAT(x) std::cout << #x << ": " << std::endl << x << std::endl
#define blockVector3(a) block<3, 1>(3 * (a), 0)
#define blockVector4(a) block<4, 1>(4 * (a), 0)

typedef double Scalar;
typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> VectorX;
typedef Eigen::SparseMatrix<Scalar> SparseMatrix;
typedef Eigen::Quaternion<Scalar> Quaternion;
typedef Eigen::Matrix<Scalar, 3, 3> Matrix3;
typedef Eigen::Matrix<Scalar, 3, 1> Vector3;
typedef Eigen::Matrix<Scalar, 4, 1> Vector4;
typedef Eigen::Triplet<Scalar> Triplet;
typedef Eigen::LLT<Eigen::Matrix<Scalar, -1, -1>> LLT;
typedef Eigen::Matrix<Scalar, 4, 4> Matrix4;

#endif //__HELPERS_PD__
