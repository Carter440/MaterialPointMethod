#include <iostream>
#include <fstream>
#include "Grain.h"
#include "MP_Grid.h"

#define TIMESTEPS 5000
#define WRITENSTEPS 10
#define DELTA_TIME 0.001
#define NUMGRAINS 2
#define GRAINSIZE 4
#define POINTSIZE 0.0625
//NOTE BULK = LAMB + (2/3)MEW
//NOTE POISSON = LAMB/(2*(LAMB+MEW)) NOTENOTE for POSSONS = 0.25 => LAMB=MEW
#define MEW 1000
#define LAMB 1000
#define EPSILON 0.0000000001
#define GRAV -10.0
#define SIZEX 22
#define SIZEY 22

using namespace std;

Matrix2d solve_stress(Matrix2d & deform, double mew, double lamb)
{//"The standard equations for stress" --some dweeb --c. schroeder --carter S
	double det = deform.determinant();
	return 4.0 * mew*(deform*deform.transpose()*deform - deform) + lamb * (det - 1.0) * det * deform.inverse().transpose();
}
//For Debug purposes only!
void solve_energy(Matrix2d & deform, double & gamma)
{
	Matrix2d I;
	I << 1.0, 0.0
		,0.0, 1.0;
	double det = deform.determinant();
	gamma = POINTSIZE * (MEW*(deform.transpose()*deform - I).cwiseAbs2().sum() + (LAMB / 2.0)*((det - 1)*(det - 1)));
}

void solve_internal_force(MP_Grid & grid, mat_point * mp)
{
	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			grid.at(mp->cell[0] + j, mp->cell[1] + i)->internal_forces -= POINTSIZE * solve_stress(mp->deformation, MEW, LAMB)*mp->deformation.transpose() * mp->interpd[(i * 2 + j)];
		}
	}
}

void update_deformation(MP_Grid & grid, mat_point * mp, double delta_t)
{
	Matrix2d I;
	I << 1.0, 0.0
		,0.0, 1.0;
	Matrix2d grad_v;
	grad_v << 0.0, 0.0
			 ,0.0, 0.0;

	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			grad_v += grid.at(mp->cell[0] + j,mp->cell[1] + i)->velocity * mp->interpd[(i * 2 + j)].transpose();
		}
	}
	mp->deformation = (I + delta_t*grad_v) * mp->deformation;

}

void apply_accelleration(MP_Grid & grid, Vector2d accell)
{
	Vector2i dim = grid.get_dimensions();
	for (int i = 0; i < dim[0]; i++)
	{
		for (int j = 0; j < dim[1]; j++)
		{
			grid.at(i, j)->external_forces += accell*grid.at(i,j)->mass;
		}
	}
}

void forward_solve(MP_Grid & grid, double delt)
{
	Vector2i dim = grid.get_dimensions();
	for (int i = 0; i < dim[0]; i++)
	{
		for (int j = 0; j < dim[1]; j++)
		{
			if (grid.at(i, j)->mass > 0.0000000001)
			{
				Vector2d accel = (grid.at(i, j)->external_forces + grid.at(i, j)->internal_forces) / grid.at(i, j)->mass;
				grid.at(i, j)->velocity = grid.at(i, j)->velocity + delt * accel;
				grid.at(i, j)->position = grid.at(i, j)->position + delt * grid.at(i, j)->velocity;
			}
		}
	}
}

void apply_boundaries(MP_Grid & grid)
{
	Vector2i dim = grid.get_dimensions();
	for (int i = 0; i < dim[0]; i++)
	{
		for (int j = 0; j < dim[1]; j++)
		{
			if (i < 2 || i > dim[0] - 2 || j < 2 || j > dim[1] - 2)
			{//Apply Boundary Conditions (FREEZE)
				grid.at(i, j)->velocity << 0.0, 0.0;
				grid.at(i,j)->position << (i - (dim[0] / 2.0)), (j - (dim[1] / 2.0));
			}
		}
	}
}

int main(int argc, char * argv[])
{

	ofstream outfile;
	outfile.open("OUTPUT.txt", std::ios_base::out);
	//should start file with num_steps num_particles
	Grain * grains[NUMGRAINS];
	int total_points = 0;
	
	Vector2d v0;
	v0 << 0.0, 5.0;
	grains[0] = new Grain(v0, 100.0, GRAINSIZE);
	for (int j = 0; j < grains[0]->get_num_points(); j++)
	{
		grains[0]->get_point(j)->velocity << 0.0, 0.0;
	}
	Vector2d v1;
	v1 << 1.0, -5.0;
	grains[1] = new Grain(v1, 100.0, GRAINSIZE);
	for (int j = 0; j < grains[1]->get_num_points(); j++)
	{
		grains[1]->get_point(j)->velocity << 0.0, 0.0;
	}
	total_points += grains[0]->get_num_points();
	total_points += grains[1]->get_num_points();
	
	MP_Grid sim_grid(SIZEX,SIZEY);

	outfile << TIMESTEPS/WRITENSTEPS << "\n";
	outfile <<  total_points << "\n";

	for (int i = 0; i < TIMESTEPS; i++)
	{
		//Extrapolate to Grid
		for (int j = 0; j < NUMGRAINS; j++)
		{
			for (int k = 0; k < grains[j]->get_num_points(); k++)
			{
				sim_grid.update_grid_position(grains[j]->get_point(k));
				sim_grid.extrapolate(grains[j]->get_point(k));
			}
		}

		sim_grid.update_velocity();

		//Update Internal Force
		for (int j = 0; j < NUMGRAINS; j++)
		{
			for (int k = 0; k < grains[j]->get_num_points(); k++)
			{
				solve_internal_force(sim_grid, grains[j]->get_point(k));
			}
		}

		//Apply external acceleration (gravity)
		apply_accelleration(sim_grid, { 0.0, GRAV });

		//Forward Euler Solve
		forward_solve(sim_grid, DELTA_TIME);

		apply_boundaries(sim_grid);

		for (int j = 0; j < NUMGRAINS; j++)
		{
			for (int k = 0; k < grains[j]->get_num_points(); k++)
			{
				//Interpolate to Points
				sim_grid.interpolate(grains[j]->get_point(k), DELTA_TIME);

				//update deformaation at point
				update_deformation(sim_grid, grains[j]->get_point(k), DELTA_TIME);
			}
		}

		//reset griid for next run
		sim_grid.reset();

		//output to file/console
		if (i % WRITENSTEPS == 0)
		{
			int point_no = 0;
			cout << "Timestep#" << i << "\n";
			outfile << "Timestep#" << i << "\n";
			for (int j = 0; j < NUMGRAINS; j++)
			{
				for (int k = 0; k < grains[j]->get_num_points(); k++)
				{
					outfile << "Point#" << point_no << "\n"
						<< grains[j]->get_point(k)->position[0] << "\n"
						<< grains[j]->get_point(k)->position[1] << "\n\n";
					point_no += 1;
				}
			}
		}
	}
	outfile.close();
	int x;
	std::cin >> x;
	return x;
}