/* R3grid.cpp
 * R3 under discretization (discretize functor) to a grid
 * Ernest Yeung  ernestyalumni@gmail.com
 * 20160630
 */
#include "R3grid.h"


Grid3d :: Grid3d(std::array<int,3> Ld_in, std::array<float,3> ld_in)
	: Ld(Ld_in), ld(ld_in)
{
	hd = { ld[0]/Ld[0], ld[1]/Ld[1], ld[2]/Ld[2]  };
	
	temperature = new float[ this->NFLAT() ];

}

std::array<float,3> Grid3d :: gridpt_to_space(std::array<int,3> index) {
	std::array<float,3> Xi { index[0]*hd[0], index[1]*hd[1], index[2]*hd[2]  } ;
	return Xi;
}

int Grid3d :: NFLAT() {
	return Ld[0]*Ld[1]*Ld[2] ;

}	

int Grid3d :: flatten(const int i_x, const int i_y, const int i_z ) {
	return i_x+i_y*Ld[0]+i_z*Ld[0]*Ld[1]  ;
}


Grid3d::~Grid3d() {
	delete[] temperature;
}

float gaussian3d(float A, float c, std::array<float,3> x_0, std::array<float,3> x) {
	return A * expf( (-1.f)*
		( (x[0]-x_0[0])*(x[0]-x_0[0]) + (x[1]-x_0[1])*(x[1]-x_0[1]) +(x[2]-x_0[2])*(x[2]-x_0[2])   ) /  
			(2.f * c*c) );
}

