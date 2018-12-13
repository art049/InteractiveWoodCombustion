/* R3grid.h
 * R3 under discretization (discretize functor) to a grid
 * Ernest Yeung  ernestyalumni@gmail.com
 * 20160701
 */
#ifndef __R3GRID_H__
#define __R3GRID_H__

#include <array> // std::array
#include <cmath> // expf

class Grid3d
{
	public : 
		std::array<int,3> Ld;
		std::array<float,3> ld;
		std::array<float,3> hd;
		
		float* temperature;
		
		Grid3d(std::array<int,3> Ld_in, std::array<float,3> ld_in);
		
		std::array<float,3> gridpt_to_space(std::array<int,3> );
		
		int NFLAT();
		
		int flatten(const int i_x, const int i_y, const int i_z ) ;
	
		~Grid3d();	
};

float gaussian3d(float A, float c, std::array<float,3> x_0, std::array<float,3> x);

#endif // __R3GRID_H__
