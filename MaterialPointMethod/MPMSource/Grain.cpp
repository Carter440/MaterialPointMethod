#include "Grain.h"
#include "Eigen/Dense"
#include <math.h>
using namespace Eigen;


Grain::Grain(Vector2d & center, double grain_mass, double radius)
{
	assert(radius > 1.0);
	int x = floor(radius);
	x = 4 * x;
	x = x * x;
	material_points =  new mat_point * [x];
	num_points = 0;
	Vector2d Lcorner;
	Lcorner << center[0]-radius, center[1]-radius;
	for (double i = 0.25; i < (2 * radius) - 0.5; i += 0.5)
	{
		for (double j = 0.25; j < (2 * radius) - 0.5; j += 0.5)
		{
			Vector2d candidate;
			candidate << Lcorner[0] + j, Lcorner[1] + i;
			//check if candidate is within grain
			if ( (candidate - center).norm() <= radius - 0.25)
			{
				material_points[num_points] = new mat_point;
				material_points[num_points]->velocity << 0.0, 0.0;
				material_points[num_points]->deformation << 1.0, 0.0, 0.0, 1.0;
				material_points[num_points]->cell << 0, 0;
				material_points[num_points]->position = candidate;
				num_points += 1;
			}
		}
	}

	for ( int q = 0; q < num_points; q++)
	{
		material_points[q]->mass = grain_mass / num_points;
	}

}

Grain::~Grain()
{
}

int Grain::get_num_points()
{
	return num_points;
}

mat_point * Grain::get_point(int index)
{
	if (index < num_points)
	{
		return material_points[index];
	}
	return nullptr;
}
