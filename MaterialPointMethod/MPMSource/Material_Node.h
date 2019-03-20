#ifndef MATERIAL_NODE

#define MATERIAL_NODE

#include <Eigen/Dense>

using namespace Eigen;

struct mat_node
{
	double mass = 0.0;
	Vector2d position;
	Vector2d velocity;
	Vector2d momentum;
	Vector2d internal_forces;
	Vector2d external_forces;
};

#endif