/* heat_3d.cu
 * 3-dim. Laplace eq. (heat eq.) by finite difference with shared memory
 * Ernest Yeung  ernestyalumni@gmail.com
 * 20160729
 */
#include "pressure.cuh"
#include "advection.cuh"
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
    if(col >= GRID_COUNT || row >= GRID_COUNT || z >= GRID_COUNT || col < 0 || row < 0 || z < 0)
        printf("flatten oob");
    return idxClip(col, GRID_COUNT) + idxClip(row,GRID_COUNT)*GRID_COUNT + idxClip(z,GRID_COUNT)*GRID_COUNT*GRID_COUNT;
}
__device__ int vflatten(int col, int row, int z) {
    if(col >= GRID_COUNT+1 || row >= GRID_COUNT+1 || z >= GRID_COUNT+1 || col < 0 || row < 0 || z < 0)
        printf("vflatten oob");
    return idxClip(col, GRID_COUNT+1) + idxClip(row,GRID_COUNT+1)*(GRID_COUNT+1) + idxClip(z,GRID_COUNT+1)*(GRID_COUNT+1)*(GRID_COUNT+1);
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
__global__ void resetKernelCentered(float* d_temp, float*  d_oldtemp, float*  d_smokedensity,float*  d_oldsmokedensity) {
    const int k_x = blockIdx.x*blockDim.x + threadIdx.x;
    const int k_y = blockIdx.y*blockDim.y + threadIdx.y;
    const int k_z = blockIdx.z*blockDim.z + threadIdx.z;
    if ((k_x >= GRID_COUNT) || (k_y >= GRID_COUNT) || (k_z >= GRID_COUNT)) return;
    
    const int k = flatten(k_x, k_y, k_z);
    d_temp[k] = d_oldtemp[k] = T_AMBIANT;
    d_smokedensity[k] = d_oldsmokedensity[k] = 0.f;
    if(d_abs(k_z - GRID_COUNT/2) * d_abs(k_z - GRID_COUNT/2) + 
    d_abs(k_y - GRID_COUNT/2) * d_abs(k_y - GRID_COUNT/2) +
    d_abs(k_x - GRID_COUNT/2) * d_abs(k_x - GRID_COUNT/2) < GRID_COUNT  *GRID_COUNT / (5*5*25)){
        d_oldsmokedensity[k] = d_smokedensity[k] =1.f;
        d_oldtemp[k] = d_temp[k] = T_AMBIANT + 50.f;
    }
}
__global__ void resetKernelVelocity(float3 *  d_vel,float3 *  d_oldvel) {
    const int k_x = blockIdx.x*blockDim.x + threadIdx.x;
    const int k_y = blockIdx.y*blockDim.y + threadIdx.y;
    const int k_z = blockIdx.z*blockDim.z + threadIdx.z;
    if ((k_x >= GRID_COUNT+1) || (k_y >= GRID_COUNT+1) || (k_z >= GRID_COUNT+1)) return;
    
    const int k = vflatten(k_x, k_y, k_z);
    d_vel[k] = make_float3(0.f, 0.f, 0.f);
    d_oldvel[k] = make_float3(0.f, 0.f, 0.f);
}


__device__ float3 getAlpham (float3 * d_vel, float3 pos, int k){
    // Iteratively compute alpha_m
    float3 alpha_m = d_vel[k] * dev_Deltat[0];
    for(uint i = 0; i < SEMILAGRANGIAN_ITERS; i++){
        float3 estimated = pos - alpha_m;
        if(estimated.x < BLOCK_SIZE) estimated.x = BLOCK_SIZE;
        if(estimated.y < BLOCK_SIZE) estimated.y = BLOCK_SIZE;
        if(estimated.z < BLOCK_SIZE) estimated.z = BLOCK_SIZE;
        if(estimated.x > GRID_SIZE - BLOCK_SIZE) estimated.x = GRID_SIZE - BLOCK_SIZE;
        if(estimated.y > GRID_SIZE - BLOCK_SIZE) estimated.y = GRID_SIZE - BLOCK_SIZE;
        if(estimated.z > GRID_SIZE - BLOCK_SIZE) estimated.z = GRID_SIZE - BLOCK_SIZE;
        uint3 b = {static_cast<uint>(estimated.x/BLOCK_SIZE),
                   static_cast<uint>(estimated.y/BLOCK_SIZE),
                   static_cast<uint>(estimated.z/BLOCK_SIZE)};
        float3 localCoord = (estimated - make_float3(b.x*BLOCK_SIZE, b.y*BLOCK_SIZE, b.z*BLOCK_SIZE)) * (1/BLOCK_SIZE); 
        alpha_m.x = (1-localCoord.x) * d_vel[vflatten(b.x, b.y, b.z)  ].x+
                    (localCoord.x)     * d_vel[vflatten(b.x+1, b.y, b.z)].x;
        alpha_m.y = (1-localCoord.y) * d_vel[vflatten(b.x, b.y, b.z)  ].y+
                    (localCoord.y)     * d_vel[vflatten(b.x, b.y+1, b.z)].y;
        alpha_m.z = (1-localCoord.z) * d_vel[vflatten(b.x, b.y, b.z)  ].z+
                    (localCoord.z)     * d_vel[vflatten(b.x, b.y, b.z+1)].z;
        alpha_m = alpha_m * dev_Deltat[0];
    }
    //CLIPPING ON FACES
    return alpha_m;
}
__device__ float3 fbuoyancy(float * d_smoke, float* d_temp, int k_x, int k_y, int k_z){
    const int k = flatten(k_x, k_y, k_z);
    float3 f = make_float3(0,0,0);
    if(k_y == GRID_COUNT - 1){
        f.y += -BUOY_ALPHA*d_smoke[k];
        f.y += BUOY_BETA*(d_temp[k] - T_AMBIANT);
    }
    else {
        f.y += -0.5 * BUOY_ALPHA*(d_smoke[k]+d_smoke[flatten(k_x,k_y+1,k_z)]);
        f.y += BUOY_BETA*((d_temp[k]+d_temp[flatten(k_x,k_y+1,k_z)]) * 0.5f - T_AMBIANT);
    }
    return f;
}
__device__ float3 fconfinement(float3 * d_vorticity, int k_x, int k_y, int k_z){
    const int k = flatten(k_x, k_y, k_z);
    if(k_x == GRID_COUNT -1 || k_y == GRID_COUNT -1 || k_z == GRID_COUNT -1)
        return make_float3(0,0,0);
    vec3 N(vec3(d_vorticity[flatten(k_x+1, k_y, k_z)]).length() - vec3(d_vorticity[k]).length(),
           vec3(d_vorticity[flatten(k_x, k_y+1, k_z)]).length() - vec3(d_vorticity[k]).length(),
           vec3(d_vorticity[flatten(k_x, k_y, k_z+1)]).length() - vec3(d_vorticity[k]).length());
    if(N.length() < 1e-6) return make_float3(0,0,0);
    N.make_unit_vector();
    vec3 f = VORTICITY_EPSILON * BLOCK_SIZE * cross(N, vec3(d_vorticity[k]));
    return f.toFloat3();
}
__global__ void computeVorticity(float3 *d_vorticity, float3* d_vel, float3* d_ccvel){
    const int k_x = threadIdx.x + blockDim.x * blockIdx.x;
    const int k_y = threadIdx.y + blockDim.y * blockIdx.y;
    const int k_z = threadIdx.z + blockDim.z * blockIdx.z;
    if ((k_x >= dev_Ld[0] ) || (k_y >= dev_Ld[1] ) || (k_z >= dev_Ld[2])) return;
    const int k = flatten(k_x, k_y, k_z);
    if(!k_x || k_x == GRID_COUNT - 1 || !k_y || k_y == GRID_COUNT -1 || !k_z || k_z == GRID_COUNT -1){
        d_vorticity[k] = make_float3(0,0,0);
        return;
    }
    d_ccvel[k].x = (d_vel[vflatten(k_x, k_y, k_z)].x + d_vel[vflatten(k_x+1, k_y, k_z)].x) * 0.5f;
    d_ccvel[k].y = (d_vel[vflatten(k_x, k_y, k_z)].y + d_vel[vflatten(k_x, k_y+1, k_z)].y) * 0.5f;
    d_ccvel[k].x = (d_vel[vflatten(k_x, k_y, k_z)].z + d_vel[vflatten(k_x, k_y, k_z+1)].z) * 0.5f;
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
__global__ void sourcesKernel(float* d_smokedensity, float* d_temp){
    const int k_x = threadIdx.x + blockDim.x * blockIdx.x;
    const int k_y = threadIdx.y + blockDim.y * blockIdx.y;
    const int k_z = threadIdx.z + blockDim.z * blockIdx.z;
    if ((k_x >= GRID_COUNT ) || (k_y >= GRID_COUNT ) || (k_z >= GRID_COUNT)) return;
    const int k = flatten(k_x, k_y, k_z);
    if(d_abs(k_z - GRID_COUNT/2) * d_abs(k_z - GRID_COUNT/2) + 
    d_abs(k_y - GRID_COUNT/2) * d_abs(k_y - GRID_COUNT/2) +
    d_abs(k_x - GRID_COUNT/2) * d_abs(k_x - GRID_COUNT/2) < GRID_COUNT  *GRID_COUNT / (5*5*25)){
        d_smokedensity[k] =1;
        d_temp[k] = T_AMBIANT + 250.f;
    }   
}

__global__ void velocityKernel(float *d_temp, float3* d_vel, float3* d_oldvel, float* d_smokedensity, float3* d_vorticity){
    const int k_x = threadIdx.x + blockDim.x * blockIdx.x;
    const int k_y = threadIdx.y + blockDim.y * blockIdx.y;
    const int k_z = threadIdx.z + blockDim.z * blockIdx.z;
    if ((k_x >= GRID_COUNT ) || (k_y >= GRID_COUNT ) || (k_z >= GRID_COUNT)) return;
    const int k = flatten(k_x, k_y, k_z);
    
    // External forces
    float3 f = {0, 0, 0};
    float3 fext = {-0.,0,0. };
    f = f + fext;
    f = f + fconfinement(d_vorticity, k_x, k_y, k_z);    
    f = f  + fbuoyancy(d_smokedensity, d_temp, k_x, k_y, k_z);
    
    //Boundary conditions
    if(k_x == 0 || k_x == GRID_COUNT - 1)
        d_vel[k].x = 0;
    if(k_y == 0 || k_y == GRID_COUNT - 1)
        d_vel[k].y = 0;
    if(k_z == 0 || k_z == GRID_COUNT - 1)
        d_vel[k].z= 0;
    
    // Semi Lagrangian Advection
    float3 pos = make_float3((k_x+0.5f)*BLOCK_SIZE, (k_y+0.5f)*BLOCK_SIZE, (k_z+0.5f)*BLOCK_SIZE);
    float3 alpha_m = getAlpham(d_oldvel, pos, k);
    // Backtracing 
    float3 estimated = pos - 2 * alpha_m;
    //Clip on boundaries faces
    if(estimated.x < BLOCK_SIZE) estimated.x = BLOCK_SIZE;
    if(estimated.y < BLOCK_SIZE) estimated.y = BLOCK_SIZE;
    if(estimated.z < BLOCK_SIZE) estimated.z = BLOCK_SIZE;
    if(estimated.x > GRID_SIZE - BLOCK_SIZE) estimated.x = GRID_SIZE - BLOCK_SIZE;
    if(estimated.y > GRID_SIZE - BLOCK_SIZE) estimated.y = GRID_SIZE - BLOCK_SIZE;
    if(estimated.z > GRID_SIZE - BLOCK_SIZE) estimated.z = GRID_SIZE - BLOCK_SIZE;
    
    int3 b = {static_cast<int>(estimated.x/BLOCK_SIZE),
               static_cast<int>(estimated.y/BLOCK_SIZE),
               static_cast<int>(estimated.z/BLOCK_SIZE)};
    
    float3 localCoord = (estimated - make_float3(b.x*BLOCK_SIZE, b.y*BLOCK_SIZE, b.z*BLOCK_SIZE)) * (1/BLOCK_SIZE);
    //Velocity per component
    float3 dv;
    dv.x = (1-localCoord.x) * d_oldvel[vflatten(b.x, b.y, b.z)  ].x+
           (localCoord.x)     * d_oldvel[vflatten(b.x+1, b.y, b.z)].x;
    dv.y = (1-localCoord.y) * d_oldvel[vflatten(b.x, b.y, b.z)  ].y+
           (localCoord.y)     * d_oldvel[vflatten(b.x, b.y+1, b.z)].y;
    dv.z = (1-localCoord.z) * d_oldvel[vflatten(b.x, b.y, b.z)  ].z+
           (localCoord.z)     * d_oldvel[vflatten(b.x, b.y, b.z+1)].z;
    //dv = dv * 2 * dev_Deltat[0];
    d_vel[k] = dv + f * dev_Deltat[0];        

    //d_vel[k] = d_oldvel[k] + dv + f * dev_Deltat[0];        
}
__device__ float scalarLinearInt(float* scalarField, float3 pos, float oobvalue){
    // Trilinear interpolation
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
    //float3 pos = make_float3((k_x)*BLOCK_SIZE, (k_y)*BLOCK_SIZE, (k_z)*BLOCK_SIZE);
    float3 alpha_m = getAlpham(d_vel, pos, k);
    // Backtracing 
    float3 estimated = pos - 2 * alpha_m;
     //Clip on boundaries faces
     if(estimated.x < BLOCK_SIZE) estimated.x = BLOCK_SIZE;
     if(estimated.y < BLOCK_SIZE) estimated.y = BLOCK_SIZE;
     if(estimated.z < BLOCK_SIZE) estimated.z = BLOCK_SIZE;
     if(estimated.x > GRID_SIZE - BLOCK_SIZE) estimated.x = GRID_SIZE - BLOCK_SIZE;
     if(estimated.y > GRID_SIZE - BLOCK_SIZE) estimated.y = GRID_SIZE - BLOCK_SIZE;
     if(estimated.z > GRID_SIZE - BLOCK_SIZE) estimated.z = GRID_SIZE - BLOCK_SIZE;
    //uint3 b = {static_cast<uint>(estimated.x/BLOCK_SIZE),
    //           static_cast<uint>(estimated.y/BLOCK_SIZE),
    //           static_cast<uint>(estimated.z/BLOCK_SIZE)};
    //float3 localCoord = (estimated - make_float3(b.x*BLOCK_SIZE, b.y*BLOCK_SIZE, b.z*BLOCK_SIZE)) * (1 / BLOCK_SIZE);
    //float ds = d_smoke[flatten(b.x, b.y, b.z) ];

    float ds = scalarLinearInt(d_oldsmoke, estimated, 0.f);
     
    
    //ds = ds * 2 * dev_Deltat[0];
    __syncthreads();
    d_smoke[k] = ds;

}
__global__ void tempAdvectionKernel(float *d_temp, float * d_oldtemp, float3* d_vel){
    const int k_x = threadIdx.x + blockDim.x * blockIdx.x;
    const int k_y = threadIdx.y + blockDim.y * blockIdx.y;
    const int k_z = threadIdx.z + blockDim.z * blockIdx.z;
    if ((k_x >= dev_Ld[0] ) || (k_y >= dev_Ld[1] ) || (k_z >= dev_Ld[2])) return;
    const int k = flatten(k_x, k_y, k_z, dev_Ld[0], dev_Ld[1],dev_Ld[2]);
    // Advection
    //float3 pos = make_float3((k_x+0.5f)*BLOCK_SIZE, (k_y+0.5f)*BLOCK_SIZE, (k_z+0.5f)*BLOCK_SIZE);
    float3 pos = make_float3((k_x)*BLOCK_SIZE, (k_y)*BLOCK_SIZE, (k_z)*BLOCK_SIZE);
    float3 alpha_m = getAlpham(d_vel, pos, k);
    // Backtracing 
    float3 estimated = pos - 2 * alpha_m;
    // Clipping
    if(estimated.x < BLOCK_SIZE) estimated.x = BLOCK_SIZE;
    if(estimated.y < BLOCK_SIZE) estimated.y = BLOCK_SIZE;
    if(estimated.z < BLOCK_SIZE) estimated.z = BLOCK_SIZE;
    if(estimated.x > GRID_SIZE - BLOCK_SIZE) estimated.x = GRID_SIZE - BLOCK_SIZE;
    if(estimated.y > GRID_SIZE - BLOCK_SIZE) estimated.y = GRID_SIZE - BLOCK_SIZE;
    if(estimated.z > GRID_SIZE - BLOCK_SIZE) estimated.z = GRID_SIZE - BLOCK_SIZE;
    //uint3 b = {static_cast<uint>(estimated.x/BLOCK_SIZE),
    //           static_cast<uint>(estimated.y/BLOCK_SIZE),
    //           static_cast<uint>(estimated.z/BLOCK_SIZE)};
    //float3 localCoord = (estimated - make_float3(b.x*BLOCK_SIZE, b.y*BLOCK_SIZE, b.z*BLOCK_SIZE)) * (1 / BLOCK_SIZE);
    //float ds = d_smoke[flatten(b.x, b.y, b.z) ];0

    float dt = scalarLinearInt(d_oldtemp, estimated, T_AMBIANT);
    //dt = dt * 2 * dev_Deltat[0];
    estimated = pos - alpha_m;
    //Clip on boundaries faces
    if(estimated.x < BLOCK_SIZE) estimated.x = BLOCK_SIZE;
    if(estimated.y < BLOCK_SIZE) estimated.y = BLOCK_SIZE;
    if(estimated.z < BLOCK_SIZE) estimated.z = BLOCK_SIZE;
    if(estimated.x > GRID_SIZE - BLOCK_SIZE) estimated.x = GRID_SIZE - BLOCK_SIZE;
    if(estimated.y > GRID_SIZE - BLOCK_SIZE) estimated.y = GRID_SIZE - BLOCK_SIZE;
    if(estimated.z > GRID_SIZE - BLOCK_SIZE) estimated.z = GRID_SIZE - BLOCK_SIZE;
    float dtR = TEMPERATURE_GAMMA * powf(scalarLinearInt(d_oldtemp, estimated, T_AMBIANT) - T_AMBIANT, 4);

    __syncthreads();
    //d_temp[k] = dt;
    d_temp[k] = dt + dtR * 2*dev_Deltat[0];
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

    // CFD
    
    
    computeVorticity<<<gridSize, M_in>>>(d_vorticity, d_oldvel, d_ccvel); 
    HANDLE_ERROR(cudaPeekAtLastError()); HANDLE_ERROR(cudaDeviceSynchronize());
    velocityKernel<<<gridSize, M_in>>>(d_oldtemp, d_vel, d_oldvel, d_oldsmokedensity, d_vorticity);
    HANDLE_ERROR(cudaPeekAtLastError()); HANDLE_ERROR(cudaDeviceSynchronize());
   
    //forceIncompressibility(d_vel, d_pressure);


    //advect(d_vel, d_oldvel, d_oldvel);
    //HANDLE_ERROR(cudaPeekAtLastError());
    //HANDLE_ERROR(cudaDeviceSynchronize());
    //advect(d_temp, d_oldtemp, d_vel, T_AMBIANT);
    //advect(d_smokedensity, d_oldsmokedensity, d_vel, 0.f);



    tempAdvectionKernel<<<gridSize, M_in>>>(d_temp, d_oldtemp, d_vel);
    HANDLE_ERROR(cudaPeekAtLastError());
    HANDLE_ERROR(cudaDeviceSynchronize());

    smokeAdvectionKernel<<<gridSize, M_in>>>(d_oldtemp, d_vel, d_smokedensity, d_oldsmokedensity);
    HANDLE_ERROR(cudaPeekAtLastError());
    HANDLE_ERROR(cudaDeviceSynchronize());

    sourcesKernel<<<gridSize, M_in>>>(d_smokedensity, d_temp);
    HANDLE_ERROR(cudaPeekAtLastError()); HANDLE_ERROR(cudaDeviceSynchronize());
    
    smokeRender(gridSize, d_out, d_smokedensity, d_smokeRadiance);
    HANDLE_ERROR(cudaPeekAtLastError());
    HANDLE_ERROR(cudaDeviceSynchronize());

    //tempKernel<<<gridSize, M_in, smSz>>>(d_temp, bc);
}

void resetVariables(float* d_temp,
                    float* d_oldtemp,
                    float3* d_vel, 
                    float3* d_oldvel, 
                    float* d_smokedensity,
                    float* d_oldsmokedensity,
                    float* d_pressure, 
                    dim3 Ld, BC bc, dim3 M_in) {
    const dim3 gridSizeC(blocksNeeded(GRID_COUNT, M_in.x), 
                         blocksNeeded(GRID_COUNT, M_in.y), 
                         blocksNeeded(GRID_COUNT, M_in.z));
    const dim3 gridSizeV(blocksNeeded(GRID_COUNT+1, M_in.x), 
                         blocksNeeded(GRID_COUNT+1, M_in.y), 
                         blocksNeeded(GRID_COUNT+1, M_in.z));
    resetKernelCentered<<<gridSizeC, M_in>>>(d_temp, d_oldtemp, d_smokedensity, d_oldsmokedensity);
    HANDLE_ERROR(cudaPeekAtLastError()); HANDLE_ERROR(cudaDeviceSynchronize());
    resetKernelVelocity<<<gridSizeV, M_in>>>(d_vel, d_oldvel);
    HANDLE_ERROR(cudaPeekAtLastError()); HANDLE_ERROR(cudaDeviceSynchronize());
    resetPressure<<<gridSizeC, M_in>>>(d_pressure);
    HANDLE_ERROR(cudaPeekAtLastError()); HANDLE_ERROR(cudaDeviceSynchronize());
}