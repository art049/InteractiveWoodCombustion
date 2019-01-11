#ifndef PRESSURE_CUH
#define PRESSURE_CUH

#include "cgls.cuh"
#include "physics.h"
#include "heat_3d.cuh"
#include "../cuda_common/errors.h"

void forceIncompressibility(float3 * d_vel, float* d_pressure);

#endif /* PRESSURE_CUH */
