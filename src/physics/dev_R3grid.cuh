/* dev_R3grid.h
 * R2 under discretization (discretize functor) to a thread block
 * Ernest Yeung  ernestyalumni@gmail.com
 * 20160728
 */
#ifndef __DEV_R3GRID_H__
#define __DEV_R3GRID_H__

#include "../cuda_common/errors.h"

extern __constant__ int dev_Ld[3]; // L^{d=2} = (L_x,L_y) \in \mathbb{N}^2

class dev_Grid3d
{
	public:
		dim3 Ld;
	
		float *dev_temperature;

		// constructor
		__host__ dev_Grid3d( dim3 );

//		__host__ ~dev_Grid3d();

		__host__ int NFLAT();
};

class dev_block3d
{
	public:
		unsigned int  N_is[3];

		__device__ dev_block3d(unsigned int N_x, unsigned int N_y, unsigned int N_z);

		__device__ int flatten(int i_x, int i_y, int i_z);
};

#endif // __DEV_R3GRID_H__
