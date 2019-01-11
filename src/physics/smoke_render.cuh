#ifndef SMOKE_RENDER_CUH
#define SMOKE_RENDER_CUH

#include "physics.h"
#include "vec3.cuh"

__host__ __device__ bool rayGridIntersect(const vec3 ray_orig, const vec3 ray_dir, 
                                 int3 * voxel, float * t);
void smokeRender(dim3 gridSize, uchar4* d_out, float * d_smokedensity, float * d_smokeRadiance);

#endif /* SMOKE_RENDER_CUH */
