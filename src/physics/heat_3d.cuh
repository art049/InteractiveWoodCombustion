/* heat_3d.h
 * 3-dim. Laplace eq. (heat eq.) by finite difference with shared memory
 * Ernest Yeung  ernestyalumni@gmail.com
 * 20160729
 */
#ifndef __HEAT_3D_H__
#define __HEAT_3D_H__


#include "../cuda_common/finitediff.cuh"
#include "dev_R3grid.cuh"

extern __constant__ float dev_Deltat[1]; // Deltat

extern __constant__ float dev_heat_params[2] ; // dev_heat_params[0] = \kappa, 
					// dev_heat_params[1] = c_V = heat capacity at constant volume per volume

struct uchar4;
// struct BC that contains all the boundary conditions
typedef struct {
	int x, y; // x and y location of pipe center
	float rad; // radius of pipe
	int chamfer; // chamfer
	float t_s, t_a, t_g; // temperatures in pipe, air, ground
} BC;

#include "physics.h"

void kernelLauncher(uchar4 *d_out,float *d_temp, float3* d_vel, float* d_smokedensity, float* smokeRadiance, dim3 Ld, BC bc, dim3 M_in, unsigned int slice) ; 


void resetVariables(float* d_temp, float3* d_vel, float* d_smokedensity, dim3 Ld, BC bc, dim3 M_in);

#endif // __HEAT_2D_H__
 
