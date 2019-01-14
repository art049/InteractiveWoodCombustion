#ifndef ADVECTION_CUH
#define ADVECTION_CUH

#include "physics.h"
#include "heat_3d.cuh"
void advect(float3 * phi, float3 * oldphi, float3 * vel);
void advect(float * phi, float * oldphi, float3 * vel, float boundary);


#endif /* ADVECTION_CUH */
