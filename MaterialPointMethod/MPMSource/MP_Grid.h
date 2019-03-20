#ifndef MP_GRID

#define MP_GRID

#include "Material_Point.h"
#include "Material_Node.h"
#include <math.h>

class MP_Grid
{
private:

	mat_node *** grid;

	int size_x;
	int size_y;

public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	MP_Grid(int dim_x, int dim_y);
	~MP_Grid();

	Vector2i get_dimensions();
	void update_grid_position(mat_point * mp);
	void update_velocity();
	void extrapolate(const mat_point * mp);
	void interpolate(mat_point * mp, double delta_t);
	void reset();

	mat_node * at(int x,int y);

};


#endif