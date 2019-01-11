/* heat_3d.cu
 * 3-dim. Laplace eq. (heat eq.) by finite difference with shared memory
 * Ernest Yeung  ernestyalumni@gmail.com
 * 20160729
 */
#include "pressure.cuh"
#include "smoke_render.cuh"
#include "heat_3d.cuh"
#include "physics.h"
#include "vec3.cuh"

#define RAD 1 // radius of the stencil; helps to deal with "boundary conditions" at (thread) block's ends

__constant__ float dev_Deltat[1];

__constant__ float dev_heat_params[2];



int blocksNeeded( int N_i, int M_i) { return (N_i+M_i-1)/M_i; }

__device__ unsigned char clip(int n) { return n > 255 ? 255 : (n < 0 ? 0 : n);}

__device__ int idxClip( int idx, int idxMax) {
    return idx > (idxMax - 1) ? (idxMax - 1): (idx < 0 ? 0 : idx);
}

__device__ int flatten(int col, int row, int z, int width, int height, int depth) {
    return idxClip(col, width) + idxClip(row,height)*width + idxClip(z,depth)*width*height;
}
__device__ int flatten(int col, int row, int z) {
    return idxClip(col, GRID_COUNT) + idxClip(row,GRID_COUNT)*GRID_COUNT + idxClip(z,GRID_COUNT)*GRID_COUNT*GRID_COUNT;
}
__device__ float3 operator+(const float3 &a, const float3 &b) {
    return make_float3(a.x+b.x, a.y+b.y, a.z+b.z);
  }
__device__ float3 operator-(const float3 &a, const float3 &b) {
    return make_float3(a.x-b.x, a.y-b.y, a.z-b.z);
  }
__device__ float3 operator*(const float3 &a, const float &b) {
  return make_float3(a.x*b, a.y*b, a.z*b);
}
__device__ float3 operator*(const float &b, const float3 &a) {
    return make_float3(a.x*b, a.y*b, a.z*b);
  }
__device__ int d_abs(int a) {
    return a > 0 ? a : -a;
}
__global__ void resetKernel(float * d_temp, float *  d_oldtemp,float3 *  d_vel,float3 *  d_oldvel,float *  d_smokedensity,float *  d_oldsmokedensity,  BC bc) {
    const int k_x = blockIdx.x*blockDim.x + threadIdx.x;
    const int k_y = blockIdx.y*blockDim.y + threadIdx.y;
    const int k_z = blockIdx.z*blockDim.z + threadIdx.z;

    if ((k_x >= dev_Ld[0]) || (k_y >= dev_Ld[1]) || (k_z >= dev_Ld[2])) return;
    const int k = flatten(k_x, k_y, k_z, dev_Ld[0], dev_Ld[1],dev_Ld[2]);
    d_temp[k] = d_oldtemp[k] = T_AMBIANT;
    d_vel[k] = d_oldvel[k] = {0.f, 0.f, 0.f};
    d_smokedensity[k] = d_oldsmokedensity[k] = 0.f;
    if(k_y < GRID_COUNT/6 && d_abs(k_z - GRID_COUNT/2) * d_abs(k_z - GRID_COUNT/2) + 
       d_abs(k_x - GRID_COUNT/2) * d_abs(k_x - GRID_COUNT/2) < GRID_COUNT  *GRID_COUNT / 25){
        d_smokedensity[k] = d_oldsmokedensity[k] =0.5f;
        d_temp[k] = d_oldtemp[k] = T_AMBIANT + 10.f;
    }
}


__device__ float3 getAlpham (float3 * d_vel, float3 pos, int k){
    // Iteratively compute alpha_m
    float3 alpha_m = d_vel[k] * dev_Deltat[0];
    for(uint i = 0; i < SEMILAGRANGIAN_ITERS; i++){
        float3 estimated = pos - alpha_m;
        if(estimated.x < 0) estimated.x = 0;
        if(estimated.y < 0) estimated.y = 0;
        if(estimated.z < 0) estimated.z = 0;
        uint3 b = {static_cast<uint>(estimated.x/BLOCK_SIZE),
                   static_cast<uint>(estimated.y/BLOCK_SIZE),
                   static_cast<uint>(estimated.z/BLOCK_SIZE)};
        float3 localCoord = (estimated - make_float3(b.x*BLOCK_SIZE, b.y*BLOCK_SIZE, b.z*BLOCK_SIZE)) * (1/BLOCK_SIZE); 
        alpha_m.x = (1-localCoord.x) * d_vel[flatten(b.x, b.y, b.z)  ].x+
                    localCoord.x     * d_vel[flatten(b.x+1, b.y, b.z)].x;
        alpha_m.y = (1-localCoord.y) * d_vel[flatten(b.x, b.y, b.z)  ].y+
                    localCoord.y     * d_vel[flatten(b.x, b.y+1, b.z)].y;
        alpha_m.z = (1-localCoord.z) * d_vel[flatten(b.x, b.y, b.z)  ].z+
                    localCoord.z     * d_vel[flatten(b.x, b.y, b.z+1)].z;
        alpha_m = alpha_m * dev_Deltat[0];
    }
    //CLIPPING ON FACES
    return alpha_m;
}
__device__ float3 fbuoyancy(float * d_smoke, float* d_temp, int k_x, int k_y, int k_z){
    const int k = flatten(k_x, k_y, k_z);
    float3 f = make_float3(0,0,0);
    f.y += -0.5 * BUOY_ALPHA*(d_smoke[k]+d_smoke[flatten(k_x,k_y+1,k_z)]);
    f.y += BUOY_BETA*((d_temp[k]+d_temp[flatten(k_x,k_y+1,k_z)]) * 0.5f - T_AMBIANT);
    return f;
}
__device__ float3 fconfinement(float3 * d_vorticity, int k_x, int k_y, int k_z){
    const int k = flatten(k_x, k_y, k_z);
    vec3 N(vec3(d_vorticity[flatten(k_x+1, k_y, k_z)]).length() - vec3(d_vorticity[k]).length(),
           vec3(d_vorticity[flatten(k_x, k_y+1, k_z)]).length() - vec3(d_vorticity[k]).length(),
           vec3(d_vorticity[flatten(k_x, k_y, k_z+1)]).length() - vec3(d_vorticity[k]).length());
    N /= BLOCK_SIZE; // NOT useful since we normalise
    //N.make_unit_vector();
    vec3 f = VORTICITY_EPSILON * BLOCK_SIZE * cross(N, vec3(d_vorticity[k]));
    return f.toFloat3();
}
__global__ void computeVorticity(float3 *d_vorticity, float3* d_vel, float3* d_ccvel){
    const int k_x = threadIdx.x + blockDim.x * blockIdx.x;
    const int k_y = threadIdx.y + blockDim.y * blockIdx.y;
    const int k_z = threadIdx.z + blockDim.z * blockIdx.z;
    if ((k_x >= dev_Ld[0] ) || (k_y >= dev_Ld[1] ) || (k_z >= dev_Ld[2])) return;
    
    const int k = flatten(k_x, k_y, k_z, dev_Ld[0], dev_Ld[1],dev_Ld[2]);
    d_ccvel[k].x = (d_vel[k].x + d_vel[flatten(k_x+1, k_y, k_z)].x) * 0.5f;
    d_ccvel[k].y = (d_vel[k].y + d_vel[flatten(k_x, k_y+1, k_z)].y) * 0.5f;
    d_ccvel[k].x = (d_vel[k].z + d_vel[flatten(k_x, k_y, k_z+1)].z) * 0.5f;
    __syncthreads();
    d_vorticity[k].x = d_ccvel[flatten(k_x, k_y+1, k_z)].z - d_ccvel[flatten(k_x, k_y-1, k_z)].z - 
                       d_ccvel[flatten(k_x, k_y, k_z+1)].y + d_ccvel[flatten(k_x, k_y, k_z+1)].y;
    d_vorticity[k].x /= 2 * BLOCK_SIZE;
    d_vorticity[k].y = d_ccvel[flatten(k_x, k_y, k_z+1)].x - d_ccvel[flatten(k_x, k_y, k_z-1)].x - 
                       d_ccvel[flatten(k_x+1, k_y, k_z)].z + d_ccvel[flatten(k_x-1, k_y, k_z)].z;
    d_vorticity[k].y /= 2 * BLOCK_SIZE;
    d_vorticity[k].z = d_ccvel[flatten(k_x+1, k_y, k_z)].y - d_ccvel[flatten(k_x-1, k_y, k_z)].y - 
                       d_ccvel[flatten(k_x, k_y+1, k_z)].x + d_ccvel[flatten(k_x, k_y-1, k_z)].x;
    d_vorticity[k].z /= 2 * BLOCK_SIZE;
}
#include <stdio.h>

__global__ void velocityKernel(float *d_temp, float3* d_vel, float3* d_oldvel, float* d_smokedensity, float3* d_vorticity){
    const int k_x = threadIdx.x + blockDim.x * blockIdx.x;
    const int k_y = threadIdx.y + blockDim.y * blockIdx.y;
    const int k_z = threadIdx.z + blockDim.z * blockIdx.z;
    if ((k_x >= dev_Ld[0] ) || (k_y >= dev_Ld[1] ) || (k_z >= dev_Ld[2])) return;
    const int k = flatten(k_x, k_y, k_z, dev_Ld[0], dev_Ld[1],dev_Ld[2]);
    if(k_x == 0 || k_x == GRID_COUNT - 1 || k_y == 0 || k_y == GRID_COUNT - 1 || k_z == 0 || k_z == GRID_COUNT - 1)
        return;

    // External forces
    float3 f = {0, 0, 0};
    float3 fext = {10,0,0};
    f = f + fext;
    f = f + fconfinement(d_vorticity, k_x, k_y, k_z);
    //printf("conf %f %f %f", fconf.x, fconf.y, fconf.z);
    f = f  + fbuoyancy(d_smokedensity, d_temp, k_x, k_y, k_z);
    d_vel[k] = d_oldvel[k] + f * dev_Deltat[0];
    
    // Semi Lagrangian Advection
    //float3 pos = make_float3(k_x*BLOCK_SIZE, k_y*BLOCK_SIZE, k_z*BLOCK_SIZE);
    float3 pos = make_float3((k_x+0.5f)*BLOCK_SIZE, (k_y+0.5f)*BLOCK_SIZE, (k_z+0.5f)*BLOCK_SIZE);

    float3 alpha_m = getAlpham(d_oldvel, pos, k);
    // Backtracing 
    float3 estimated = pos - 2 * alpha_m;
    if(estimated.x < 0) estimated.x = 0;
    if(estimated.y < 0) estimated.y = 0;
    if(estimated.z < 0) estimated.z = 0;
    uint3 b = {static_cast<uint>(estimated.x/BLOCK_SIZE),
               static_cast<uint>(estimated.y/BLOCK_SIZE),
               static_cast<uint>(estimated.z/BLOCK_SIZE)};
    float3 localCoord = (estimated - make_float3(b.x*BLOCK_SIZE, b.y*BLOCK_SIZE, b.z*BLOCK_SIZE)) * (1/BLOCK_SIZE);
    //Velocity per component
    float3 dv;
    dv.x = (1-localCoord.x) * d_oldvel[flatten(b.x, b.y, b.z)  ].x+
           localCoord.x     * d_oldvel[flatten(b.x+1, b.y, b.z)].x;
    dv.y = (1-localCoord.y) * d_oldvel[flatten(b.x, b.y, b.z)  ].y+
           localCoord.y     * d_oldvel[flatten(b.x, b.y+1, b.z)].y;
    dv.z = (1-localCoord.z) * d_oldvel[flatten(b.x, b.y, b.z)  ].z+
           localCoord.z     * d_oldvel[flatten(b.x, b.y, b.z+1)].z;
    dv = dv * 2 * dev_Deltat[0];
    d_vel[k] = d_vel[k] + dv;        
}
__device__ float scalarLinearInt(float* scalarField, float3 pos, float oobvalue){
    int x = static_cast<int> (pos.x / BLOCK_SIZE);
    int y = static_cast<int> (pos.y / BLOCK_SIZE);
    int z = static_cast<int> (pos.z / BLOCK_SIZE);
    //Getting voxel edges
    if(fabs(pos.x - x * BLOCK_SIZE) > fabs(pos.x - (x+1)*BLOCK_SIZE)) x++;
    if(fabs(pos.y - y * BLOCK_SIZE) > fabs(pos.y - (y+1)*BLOCK_SIZE)) y++;
    if(fabs(pos.z - z * BLOCK_SIZE) > fabs(pos.z - (z+1)*BLOCK_SIZE)) z++;
    //pos is inside voxels [x-1, x] [y-1, y] [z-1, z]
    //BOUND CHECK
    if(x <= 0 || x >= GRID_COUNT || y <= 0 || y >= GRID_COUNT || z <= 0 || z >= GRID_COUNT)
        return oobvalue;
    
    float tx = (pos.x /BLOCK_SIZE - (x - 0.5f) );
    float ty = (pos.y /BLOCK_SIZE - (y - 0.5f) );
    float tz = (pos.z /BLOCK_SIZE - (z - 0.5f) );

    // Bottom z then upper z
    float bybz = tx * scalarField[flatten(x,y-1,z-1)] + (1-tx) * scalarField[flatten(x-1,y-1,z-1)];
    float uybz = tx * scalarField[flatten(x,y,z-1)] + (1-tx) * scalarField[flatten(x-1,y,z-1)];
    float bz = (1-ty) * bybz + ty * uybz;
    float byuz = tx * scalarField[flatten(x,y-1,z)] + (1-tx) * scalarField[flatten(x-1,y-1,z)];
    float uyuz = tx * scalarField[flatten(x,y,z)] + (1-tx) * scalarField[flatten(x-1,y,z)];
    float uz = (1-ty) * byuz + ty * uyuz;
    return (1-tz) * bz + tz * uz;
}
__global__ void smokeAdvectionKernel(float *d_temp, float3* d_vel, float* d_smoke, float* d_oldsmoke){
    const int k_x = threadIdx.x + blockDim.x * blockIdx.x;
    const int k_y = threadIdx.y + blockDim.y * blockIdx.y;
    const int k_z = threadIdx.z + blockDim.z * blockIdx.z;
    if ((k_x >= dev_Ld[0] ) || (k_y >= dev_Ld[1] ) || (k_z >= dev_Ld[2])) return;
    const int k = flatten(k_x, k_y, k_z, dev_Ld[0], dev_Ld[1],dev_Ld[2]);
    // Advection
    float3 pos = make_float3((k_x+0.5f)*BLOCK_SIZE, (k_y+0.5f)*BLOCK_SIZE, (k_z+0.5f)*BLOCK_SIZE);
    float3 alpha_m = getAlpham(d_vel, pos, k);
    // Backtracing 
    float3 estimated = pos - 2 * alpha_m;
    if(estimated.x < 0) estimated.x = 0;
    if(estimated.y < 0) estimated.y = 0;
    if(estimated.z < 0) estimated.z = 0;
    //uint3 b = {static_cast<uint>(estimated.x/BLOCK_SIZE),
    //           static_cast<uint>(estimated.y/BLOCK_SIZE),
    //           static_cast<uint>(estimated.z/BLOCK_SIZE)};
    //float3 localCoord = (estimated - make_float3(b.x*BLOCK_SIZE, b.y*BLOCK_SIZE, b.z*BLOCK_SIZE)) * (1 / BLOCK_SIZE);
    //float ds = d_smoke[flatten(b.x, b.y, b.z) ];

    float ds = scalarLinearInt(d_smoke, estimated, 0.f);
    ds = ds * 2 * dev_Deltat[0];
    //NEED OLD GRID FOR THIS
    d_smoke[k] = d_oldsmoke[k] + ds;
}
__global__ void tempAdvectionKernel(float *d_temp, float * d_oldtemp, float3* d_vel, float* d_smoke){
    const int k_x = threadIdx.x + blockDim.x * blockIdx.x;
    const int k_y = threadIdx.y + blockDim.y * blockIdx.y;
    const int k_z = threadIdx.z + blockDim.z * blockIdx.z;
    if ((k_x >= dev_Ld[0] ) || (k_y >= dev_Ld[1] ) || (k_z >= dev_Ld[2])) return;
    const int k = flatten(k_x, k_y, k_z, dev_Ld[0], dev_Ld[1],dev_Ld[2]);
    // Advection
    float3 pos = make_float3((k_x+0.5f)*BLOCK_SIZE, (k_y+0.5f)*BLOCK_SIZE, (k_z+0.5f)*BLOCK_SIZE);
    float3 alpha_m = getAlpham(d_vel, pos, k);
    // Backtracing 
    float3 estimated = pos - 2 * alpha_m;
    if(estimated.x < 0) estimated.x = 0;
    if(estimated.y < 0) estimated.y = 0;
    if(estimated.z < 0) estimated.z = 0;
    //uint3 b = {static_cast<uint>(estimated.x/BLOCK_SIZE),
    //           static_cast<uint>(estimated.y/BLOCK_SIZE),
    //           static_cast<uint>(estimated.z/BLOCK_SIZE)};
    //float3 localCoord = (estimated - make_float3(b.x*BLOCK_SIZE, b.y*BLOCK_SIZE, b.z*BLOCK_SIZE)) * (1 / BLOCK_SIZE);
    //float ds = d_smoke[flatten(b.x, b.y, b.z) ];

    float dt = scalarLinearInt(d_temp, estimated, T_AMBIANT);
    dt = dt * 2 * dev_Deltat[0];
    __syncthreads();
    d_temp[k] = d_oldtemp[k] + dt;
}







void kernelLauncher(uchar4 *d_out,
                    float *d_temp, 
                    float *d_oldtemp, 
                    float3* d_vel, 
                    float3* d_oldvel, 
                    float* d_pressure,
                    float3* d_ccvel,
                    float3* d_vorticity,
                    float* d_smokedensity,
                    float* d_oldsmokedensity,
                    float * d_smokeRadiance,
                    int activeBuffer, dim3 Ld, BC bc, dim3 M_in, unsigned int slice) {
    const dim3 gridSize(blocksNeeded(Ld.x, M_in.x), blocksNeeded(Ld.y, M_in.y), 
                        blocksNeeded(Ld.z,M_in.z));
    const size_t smSz = (M_in.x + 2 * RAD)*(M_in.y + 2 * RAD)*(M_in.z + 2 * RAD)*sizeof(float);//shared mem size
    // CFD
    computeVorticity<<<gridSize, M_in>>>(d_vorticity, d_oldvel, d_ccvel);
    HANDLE_ERROR(cudaPeekAtLastError());
    velocityKernel<<<gridSize, M_in, smSz>>>(d_oldtemp, d_vel, d_oldvel, d_oldsmokedensity, d_vorticity);
    HANDLE_ERROR(cudaPeekAtLastError());
    forceIncompressibility(d_vel, d_pressure);
    tempAdvectionKernel<<<gridSize, M_in, smSz>>>(d_temp, d_oldtemp, d_oldvel, d_oldsmokedensity);
    HANDLE_ERROR(cudaPeekAtLastError());
    smokeAdvectionKernel<<<gridSize, M_in, smSz>>>(d_oldtemp, d_oldvel, d_smokedensity, d_oldsmokedensity);
    HANDLE_ERROR(cudaPeekAtLastError());
    
    smokeRender(gridSize, d_out, d_smokedensity, d_smokeRadiance);

    //tempKernel<<<gridSize, M_in, smSz>>>(d_temp, bc);
    HANDLE_ERROR(cudaDeviceSynchronize());
}

void resetVariables(float* d_temp,
                    float* d_oldtemp,
                    float3* d_vel, 
                    float3* d_oldvel, 
                    float* d_smokedensity,
                    float* d_oldsmokedensity, 
                    dim3 Ld, BC bc, dim3 M_in) {
    const dim3 gridSize( blocksNeeded(Ld.x, M_in.x), blocksNeeded( Ld.y, M_in.y), 
                            blocksNeeded(Ld.z, M_in.z));
    resetKernel<<<gridSize, M_in>>>(d_temp, d_oldtemp, d_vel, d_oldvel, d_smokedensity, d_oldsmokedensity, bc);
    HANDLE_ERROR(cudaPeekAtLastError());
    HANDLE_ERROR(cudaDeviceSynchronize());
}