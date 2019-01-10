#ifndef PRESSURE_CUH
#define PRESSURE_CUH

#include "physics.h"
#include "heat_3d.cuh"
#include "cgls.cuh"
#include "../cuda_common/errors.h"

void forceIncompressibility(float3 * d_vel, float* d_pressure);

#endif /* PRESSURE_CUH */
