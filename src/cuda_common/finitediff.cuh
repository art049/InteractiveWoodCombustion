/* finitediff.h
 * finite difference methods on a grid
 * Ernest Yeung  ernestyalumni@gmail.com
 * 20160710
 */
#ifndef __FINITEDIFF_H__
#define __FINITEDIFF_H__ 

#include "errors.h"

// I fixed the size of dev_c to be 4 because I don't anticipate that we'd need more accuracy 
extern __constant__ float3 dev_cnus[4]; // i = x,y,z; j = 1,2,3,4

void set1DerivativeParameters(const float hd_i[3] );

void set2DerivativeParameters(const float hd_i[3] );

void set3DerivativeParameters(const float hd_i[3] );

void set4DerivativeParameters(const float hd_i[3] );


__device__ float dev_dirdblder1(float centerval, float stencil[1][2], float c_nus[4]);

__device__ float dev_dirdblder2(float centerval, float stencil[2][2], float c_nus[4]);

__device__ float dev_dirdblder3(float centerval, float stencil[3][2], float c_nus[4]);

__device__ float dev_dirdblder4(float centerval, float stencil[4][2], float c_nus[4]);

// LAPLACIAN (LAP) in 1, 2, 3, 4 stencils, central difference

__device__ float dev_lap1( float centerval, float3 stencil[1][2]  ) ;

__device__ float dev_lap2( float centerval, float3 stencil[2][2]  ) ; 

__device__ float dev_lap3( float centerval, float3 stencil[3][2]  ) ;

__device__ float dev_lap4( float centerval, float3 stencil[4][2]  ) ;


#endif // __FINITEDIFF_H__
