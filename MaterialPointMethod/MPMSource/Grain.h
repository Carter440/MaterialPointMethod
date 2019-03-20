#ifndef GRAIN

#define GRAIN
#define PI 	3.14159265358979323846
#include "Material_Point.h"


class Grain
{
private:
	int num_points = 0;
	mat_point ** material_points;
	double radius;

public:
	Grain(Vector2d & center, double mass, double radius);
	~Grain();
	int get_num_points();
	mat_point * get_point(int index);

};

#endif