#ifndef MATERIAL_POINT

#define MATERIAL_POINT

#include <Eigen/Dense>

using namespace Eigen;

struct mat_point
{
	double mass = 0.0;
	Vector2d position;
	Vector2d velocity;
	double interpf[4];
	Vector2d interpd[4];
	Vector2i cell;
	Matrix2d deformation;

};

#endif