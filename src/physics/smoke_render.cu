#include "smoke_render.cuh"
#include "heat_3d.cuh"
#include "physics.h"

__host__ __device__ void swap(float * a, float * b){
    float temp = *b;
    *b = *a;
    *a = temp;
}
__host__ __device__ float posclip(float a){
    return a>0 ? a : 0;
}
__host__ __device__ bool rayGridIntersect(const vec3 ray_orig, const vec3 ray_dir, 
                                 int3 * voxel, float * t){
    const float m = 0., M = GRID_SIZE;
    float tmin = (m - ray_orig.x()) / ray_dir.x();
    float tmax = (M - ray_orig.x()) / ray_dir.x();
    if (tmin > tmax) swap(&tmin, &tmax);

    float tymin = (m - ray_orig.y()) / ray_dir.y();
    float tymax = (M - ray_orig.y()) / ray_dir.y();
    if (tymin > tymax) swap(&tymin, &tymax);
    if ((tmin > tymax) || (tymin > tmax)) return false;

    if (tymin > tmin) tmin = tymin;
    if (tymax < tmax) tmax = tymax;
    float tzmin = (m - ray_orig.z()) / ray_dir.z();
    float tzmax = (M - ray_orig.z()) / ray_dir.z();

    if (tzmin > tzmax) swap(&tzmin, &tzmax);
    if ((tmin > tzmax) || (tzmin > tmax)) return false;

    if (tzmin > tmin) tmin = tzmin;
    if (tzmax < tmax) tmax = tzmax;
    
    *t = tmin;
    if(*t < 0) return false;
    voxel->x = static_cast<int> (posclip(ray_orig.x() + *t * ray_dir.x()) / BLOCK_SIZE);
    voxel->y = static_cast<int> (posclip(ray_orig.y() + *t * ray_dir.y()) / BLOCK_SIZE);
    voxel->z = static_cast<int> (posclip(ray_orig.z() + *t * ray_dir.z()) / BLOCK_SIZE);
    return true; 

}

__global__ void smokeLightKernel(vec3 ray_o, vec3 ray_dir, float * d_smoke, float * voxelRadiance){
    const int k_x = threadIdx.x + blockDim.x * blockIdx.x;
    const int k_y = threadIdx.y + blockDim.y * blockIdx.y;
    if(k_x >= SMOKE_RAY_SQRT_COUNT || k_y >= SMOKE_RAY_SQRT_COUNT) return;

    vec3 ray_orig = ray_o + vec3(0, (k_x+0.5) * BLOCK_SIZE, (k_y+0.5) * BLOCK_SIZE);
    float ray_transparency = 1.f;
    int3 step;
    step.x = ray_dir.x() > 0 ? 1 : -1;
    step.y = ray_dir.y() > 0 ? 1 : -1;
    step.z = ray_dir.z() > 0 ? 1 : -1;
    int3 voxel; //Current voxel coordinate
    float t; // Initial hit param value
    if(!rayGridIntersect(ray_orig, ray_dir, &voxel, &t)) return;
    //printf("%d %d %d\n", voxel.x, voxel.y, voxel.z);
    float3 tMax; 
    // DIV by zero handling
    float3 tDelta = {BLOCK_SIZE / ray_dir.x(), BLOCK_SIZE / ray_dir.y(), BLOCK_SIZE / ray_dir.z()};
    while(true){
        const int k = flatten(voxel.x, voxel.y, voxel.z);
        const float voxelTransp = expf(-(SMOKE_EXTINCTION_COEFF*BLOCK_SIZE/d_smoke[k]));
        //voxelRadiance[k] += ray_transparency;
        voxelRadiance[k] += SMOKE_ALBEDO * SMOKE_LIGHT_RADIANCE * voxelTransp * ray_transparency;
        ray_transparency *= 1-voxelTransp;

        tMax.x = (BLOCK_SIZE - fmod((ray_orig + t*ray_dir).x(), BLOCK_SIZE))/ray_dir.x();
        tMax.y = (BLOCK_SIZE - fmod((ray_orig + t*ray_dir).y(), BLOCK_SIZE))/ray_dir.y();
        tMax.z = (BLOCK_SIZE - fmod((ray_orig + t*ray_dir).z(), BLOCK_SIZE))/ray_dir.z();
        if(tMax.x < tMax.y) {
            if(tMax.x < tMax.z) {
                voxel.x += step.x;
                if(voxel.x >= GRID_COUNT || voxel.x < 0)return; /* outside grid */
                tMax.x += tDelta.x;
            } else {
                voxel.z += step.z;
                if(voxel.z >= GRID_COUNT || voxel.z < 0)return;
                tMax.z += tDelta.z;
            }
        } else {
            if(tMax.y < tMax.z) {
                voxel.y += step.y;
                if(voxel.y >= GRID_COUNT || voxel.y < 0)return;
                tMax.y += tDelta.y;
            } else {
                voxel.z += step.z;
                if(voxel.z >= GRID_COUNT || voxel.z < 0)return;
                tMax.z += tDelta.z;
            }
        }
    __syncthreads();
    }
    
}
__global__ void resetSmokeRadiance(float * voxelRadiance){
    const int k_x = threadIdx.x + blockDim.x * blockIdx.x;
    const int k_y = threadIdx.y + blockDim.y * blockIdx.y;
    const int k_z = threadIdx.z + blockDim.z * blockIdx.z;
    if ((k_x >= dev_Ld[0] ) || (k_y >= dev_Ld[1] ) || (k_z >= dev_Ld[2])) return;
    const int k = flatten(k_x, k_y, k_z, dev_Ld[0], dev_Ld[1],dev_Ld[2]);
    
    voxelRadiance[k] = 0.f;
}
__global__ void generateSmokeColorBuffer( uchar4* dev_out, const float* d_smoke, float* d_smokeRadiance) {
    const int k_x = threadIdx.x + blockDim.x * blockIdx.x;
    const int k_y = threadIdx.y + blockDim.y * blockIdx.y;
    const int k_z = threadIdx.z + blockDim.z * blockIdx.z;
    if ((k_x >= dev_Ld[0] ) || (k_y >= dev_Ld[1] ) || (k_z >= dev_Ld[2])) return;
    const int k = flatten(k_x, k_y, k_z, dev_Ld[0], dev_Ld[1],dev_Ld[2]);
    if(isnan(d_smokeRadiance[k]) || isinf(d_smokeRadiance[k])) d_smokeRadiance[k] = 0;
    const unsigned char transparency = 255;// clip((int) (expf(-(fabsf(SMOKE_EXTINCTION_COEFF/d_smoke[k]))* BLOCK_SIZE)*255.f));
    const unsigned char intensity = clip((int) (d_smokeRadiance[k]*255.f));
    
    if(intensity == 0 && transparency==255)
        printf("radiance %f %d smoke %f\n", d_smokeRadiance[k], intensity,d_smoke[k]);
    //const unsigned char intensity = clip((int) (d_smoke[k]*255.f));
    for(uint i = 0; i < 4; i++){
        dev_out[4*k+i].x = intensity;
        dev_out[4*k+i].z = intensity;
        dev_out[4*k+i].y = intensity;
        dev_out[4*k+i].w = transparency; // 255 => solid display
    }
}

void smokeRender(dim3 gridSize, uchar4* d_out, float * d_smokedensity, float * d_smokeRadiance){
    // Rendering computations
    const dim3 rayBlockSize(8,8);
    const dim3 rayGridSize(blocksNeeded(SMOKE_RAY_SQRT_COUNT, rayBlockSize.x), 
                           blocksNeeded(SMOKE_RAY_SQRT_COUNT, rayBlockSize.y));
    resetSmokeRadiance<<<gridSize, M_i>>>(d_smokeRadiance);
    HANDLE_ERROR(cudaPeekAtLastError());
    HANDLE_ERROR(cudaDeviceSynchronize());
    smokeLightKernel<<<rayGridSize, rayBlockSize>>>(vec3(SMOKE_LIGHT_POS), vec3(SMOKE_LIGHT_DIR), d_smokedensity, d_smokeRadiance);
    HANDLE_ERROR(cudaPeekAtLastError());
    HANDLE_ERROR(cudaDeviceSynchronize());

    generateSmokeColorBuffer<<<gridSize, M_i>>>(d_out, d_smokedensity, d_smokeRadiance);
    HANDLE_ERROR(cudaPeekAtLastError());
    HANDLE_ERROR(cudaDeviceSynchronize());

}