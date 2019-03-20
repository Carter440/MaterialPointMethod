#include "MP_Grid.h"
#include <iostream>
#include <math.h>

MP_Grid::MP_Grid(int dim_x,int dim_y)
{
	size_x = dim_x;
	size_y = dim_y;
	grid = new __declspec(align(16)) mat_node ** [size_x];
	for (int i = 0; i < size_x; i++)
	{
		grid[i] = new mat_node * [size_y];
		for (int j = 0; j < size_y; j++)
		{
			grid[i][j] = new mat_node;
			grid[i][j]->mass = 0.0;
			grid[i][j]->velocity << 0.0, 0.0;
			grid[i][j]->momentum << 0.0, 0.0;
			grid[i][j]->internal_forces << 0.0, 0.0;
			grid[i][j]->external_forces << 0.0, 0.0;
			grid[i][j]->position << (i - (size_x/2.0)), (j - (size_y/2.0));
		}
	}

}

MP_Grid::~MP_Grid()
{

}

Vector2i MP_Grid::get_dimensions()
{
	return Vector2i(size_x, size_y);
}

void MP_Grid::update_grid_position(mat_point * mp)
{
	double gpdx;
	double gpdy;

	Vector2d pos_in_cell;
	pos_in_cell << modf(mp->position[0], &gpdx), modf(mp->position[1], &gpdy);

	if (pos_in_cell[0] < 0.0)
	{
		pos_in_cell[0] = 1.0 - abs(pos_in_cell[0]);
		gpdx -= 1.0;
	}
	if (pos_in_cell[1] < 0.0)
	{
		pos_in_cell[1] = 1.0 - abs(pos_in_cell[1]);
		gpdy -= 1.0;
	}
	mp->cell << floor(gpdx + (size_x / 2)), floor(gpdy + (size_y / 2));
	//interpolation weights
	mp->interpf[0] = (1.0 - pos_in_cell[0])*(1.0 - pos_in_cell[1]);
	mp->interpf[1] = (pos_in_cell[0])*(1.0 - pos_in_cell[1]);
	mp->interpf[2] = (1.0 - pos_in_cell[0])*(pos_in_cell[1]);
	mp->interpf[3] = (pos_in_cell[0])*(pos_in_cell[1]);

	//gradient of interpolation weights
	mp->interpd[0] << pos_in_cell[1] - 1, pos_in_cell[0] - 1;
	mp->interpd[1] << 1 - pos_in_cell[1], -pos_in_cell[0];
	mp->interpd[2] << -pos_in_cell[1], 1 - pos_in_cell[0];
	mp->interpd[3] << pos_in_cell[1], pos_in_cell[0];

}

void MP_Grid::update_velocity()
{
	for (int i = 0; i < size_x; i++)
	{
		for (int j = 0; j < size_y; j++)
		{
			if (grid[i][j]->mass > 0.0000000001)
			{
				grid[i][j]->velocity = grid[i][j]->momentum / grid[i][j]->mass;
			}
		}
	}

}

void MP_Grid::extrapolate(const mat_point * mp)
{

	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			grid[mp->cell[0] + j][mp->cell[1] + i]->mass += mp->mass * mp->interpf[(i * 2 + j)];
			grid[mp->cell[0] + j][mp->cell[1] + i]->momentum += mp->velocity * mp->mass * mp->interpf[(i * 2 + j)];
		}
	}
}
//interpolate values back to mps and clear grid point if all values interpolated back;
void MP_Grid::interpolate(mat_point * mp, double delta_t)
{

	mp->position << 0.0, 0.0;

	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			if (grid[mp->cell[0] + j][mp->cell[1] + i]->mass > 0.0000001)
			{
				mp->position += (grid[mp->cell[0] + j][mp->cell[1] + i]->position) * mp->interpf[(i * 2 + j)];
				mp->velocity += ((grid[mp->cell[0] + j][mp->cell[1] + i]->external_forces
					+ grid[mp->cell[0] + j][mp->cell[1] + i]->internal_forces)
					/ grid[mp->cell[0] + j][mp->cell[1] + i]->mass)*mp->interpf[(i * 2 + j)] * delta_t;
			}
		}
	}

}

void MP_Grid::reset()
{
	for (int i = 0; i < size_x; i++)
	{
		for (int j = 0; j < size_y; j++)
		{
			grid[i][j]->mass = 0.0;
			grid[i][j]->velocity << 0.0, 0.0;
			grid[i][j]->momentum << 0.0, 0.0;
			grid[i][j]->internal_forces << 0.0, 0.0;
			grid[i][j]->external_forces << 0.0, 0.0;
			grid[i][j]->position << (i - (size_x/2.0)), (j - (size_y/2.0));
		}
	}

}

mat_node * MP_Grid::at(int x, int y)
{
	return grid[x][y];
}
