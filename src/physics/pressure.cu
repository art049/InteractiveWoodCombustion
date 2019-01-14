#include "pressure.cuh"

__device__ int valIndex(int k_x, int k_y, int k_z){
    const int center_nnz = 7;
    const int boundary_nnz = 1;
    //INTERMEDIATE FACE SIZE
    const int cfnnz = 4*(GRID_COUNT-1)*boundary_nnz + (GRID_COUNT-2)*(GRID_COUNT-2)*center_nnz;

    if(k_z == 0){
        return (k_x + k_y*GRID_COUNT) * boundary_nnz;
    }
    else if(k_z == GRID_COUNT - 1){
        return GRID_COUNT * GRID_COUNT * boundary_nnz +
               (GRID_COUNT - 2) * cfnnz +
               (k_x + k_y*GRID_COUNT) * boundary_nnz;
    }
    else {
        if(k_y == 0){
            return GRID_COUNT * GRID_COUNT * boundary_nnz + 
                   (k_z - 1) * cfnnz + 
                   k_x * boundary_nnz;
        }
        else if(k_y == GRID_COUNT - 1){
            return GRID_COUNT * GRID_COUNT * boundary_nnz + 
                   (k_z - 1) * cfnnz + 
                   GRID_COUNT * boundary_nnz +
                   (GRID_COUNT-2) * (2*boundary_nnz + (GRID_COUNT-2)*center_nnz) +
                   (k_x) * boundary_nnz;
        }
        else{
            if(k_x == 0){
                return GRID_COUNT * GRID_COUNT * boundary_nnz +
                       (k_z - 1) * cfnnz +
                       GRID_COUNT * boundary_nnz +
                       (k_y - 1) * (2*boundary_nnz + (GRID_COUNT-2)*center_nnz);
            }
            else if(k_x == GRID_COUNT - 1){
                return GRID_COUNT * GRID_COUNT * boundary_nnz +
                       (k_z - 1) * cfnnz +
                       GRID_COUNT * boundary_nnz +
                       (k_y - 1) * (2*boundary_nnz + (GRID_COUNT-2)*center_nnz)+
                       boundary_nnz + (GRID_COUNT-2)*center_nnz;
            }
            else {
                return GRID_COUNT * GRID_COUNT * boundary_nnz +
                       (k_z - 1) * cfnnz +
                       GRID_COUNT * boundary_nnz +
                       (k_y - 1) * (2*boundary_nnz + (GRID_COUNT-2)*center_nnz)+
                       boundary_nnz + (k_x-1) * center_nnz;
            }
        }
    }
}

__global__ void prepareSystem(const int NFLAT, float3* d_vel, float* d_b, float* d_val, int* d_cind, int * d_rptr) {
    const int k_x = threadIdx.x + blockDim.x * blockIdx.x;
    const int k_y = threadIdx.y + blockDim.y * blockIdx.y;
    const int k_z = threadIdx.z + blockDim.z * blockIdx.z;
    if ((k_x >= dev_Ld[0] ) || (k_y >= dev_Ld[1] ) || (k_z >= dev_Ld[2])) return;
    const int k = flatten(k_x, k_y, k_z);
    // Matrix 
    const int offset = 7 * flatten(k_x, k_y, k_z);
    if(k_x > 0 && k_x < GRID_COUNT - 1 && k_y > 0 && k_y < GRID_COUNT - 1 && k_z > 0 && k_z < GRID_COUNT - 1) {
        //B term
        d_b[k] = d_vel[flatten(k_x+1, k_y, k_z)].x - d_vel[flatten(k_x, k_y, k_z)].x + 
                 d_vel[flatten(k_x, k_y+1, k_z)].y - d_vel[flatten(k_x, k_y, k_z)].y + 
                 d_vel[flatten(k_x, k_y, k_z+1)].z - d_vel[flatten(k_x, k_y, k_z)].z ;
        d_b[k] /= BLOCK_SIZE * dev_Deltat[0];
        d_val [offset    ] =  1;
        d_cind[offset    ] = flatten(k_x, k_y, k_z-1);
        d_val [offset + 1] =  1; 
        d_cind[offset + 1] = flatten(k_x, k_y-1, k_z);
        d_val [offset + 2] =  1; 
        d_cind[offset + 2] = flatten(k_x-1, k_y, k_z);
        d_val [offset + 3] = -6; 
        d_cind[offset + 3] = k;
        d_val [offset + 4] =  1; 
        d_cind[offset + 4] = flatten(k_x+1, k_y, k_z);
        d_val [offset + 5] =  1; 
        d_cind[offset + 5] = flatten(k_x, k_y+1, k_z);
        d_val [offset + 6] =  1; 
        d_cind[offset + 6] = flatten(k_x, k_y, k_z+1);
        
        d_rptr[k] = offset;
    }
    //PRESSURE DIRICHLET BOUNDARY CONDITION
    else {
        d_b[k] = P_ATM;
        // DUMMY VALUES to preserve alignement
        d_cind[offset    ] =  0;//(k+6)%NFLAT;
        d_cind[offset + 1] =  0;//(k+1)%NFLAT;
        d_cind[offset + 2] =  0;//(k+2)%NFLAT;
        d_cind[offset + 3] =  0;//(k+3)%NFLAT;
        d_cind[offset + 4] =  0;//(k+4)%NFLAT;
        d_cind[offset + 5] =  0;//(k+5)%NFLAT;
        d_val [offset    ] =  0;
        d_val [offset + 1] =  0; 
        d_val [offset + 2] =  0; 
        d_val [offset + 3] =  0; 
        d_val [offset + 4] =  0; 
        d_val [offset + 5] =  0; 
        d_val [offset + 6] = 1; 
        d_cind[offset + 6] = k;
        
        d_rptr[k] = offset;
        if(k == NFLAT-1){
            d_rptr[NFLAT] = 7 * NFLAT;
        }
    }
}

__global__ void prepareJacobiMethod(float3* d_vel, float* d_f) {
    const int k_x = threadIdx.x + blockDim.x * blockIdx.x;
    const int k_y = threadIdx.y + blockDim.y * blockIdx.y;
    const int k_z = threadIdx.z + blockDim.z * blockIdx.z;
    if ((k_x >= dev_Ld[0] ) || (k_y >= dev_Ld[1] ) || (k_z >= dev_Ld[2])) return;
    const int k = flatten(k_x, k_y, k_z);
    if(k_x > 0 && k_x < GRID_COUNT - 1 && k_y > 0 && k_y < GRID_COUNT - 1 && k_z > 0 && k_z < GRID_COUNT - 1){
        d_f[k] = d_vel[vflatten(k_x+1, k_y, k_z)].x - d_vel[vflatten(k_x, k_y, k_z)].x + 
                 d_vel[vflatten(k_x, k_y+1, k_z)].y - d_vel[vflatten(k_x, k_y, k_z)].y + 
                 d_vel[vflatten(k_x, k_y, k_z+1)].z - d_vel[vflatten(k_x, k_y, k_z)].z ;
        d_f[k] /= BLOCK_SIZE * dev_Deltat[0];
    }
    else
        d_f[k] = 0;
}

__global__ void jacobiIterations(float * d_pressure, float * d_temppressure, float * d_f) {
    const int k_x = threadIdx.x + blockDim.x * blockIdx.x;
    const int k_y = threadIdx.y + blockDim.y * blockIdx.y;
    const int k_z = threadIdx.z + blockDim.z * blockIdx.z;
    if(k_x == 0 || k_x >= GRID_COUNT-1  || k_y == 0 || k_y >= GRID_COUNT-1 || k_z == 0 || k_z >= GRID_COUNT-1) return;
    const int k = flatten(k_x, k_y, k_z);
    //d_pressure[k] = d_f[k];
    float * d_oldpressure = d_pressure;
    float * d_newpressure = d_temppressure;
    float * tmp;
    __syncthreads();
    for(int i=0; i < PRESSURE_JACOBI_ITERATIONS; i++){
        d_newpressure[k] = d_oldpressure[flatten(k_x, k_y, k_z-1)] +
                           d_oldpressure[flatten(k_x, k_y-1, k_z)] +
                           d_oldpressure[flatten(k_x-1, k_y, k_z)] +
                           d_oldpressure[flatten(k_x+1, k_y, k_z)] +
                           d_oldpressure[flatten(k_x, k_y+1, k_z)] +
                           d_oldpressure[flatten(k_x, k_y, k_z+1)] ;
        d_newpressure[k] /= 6;
        d_newpressure[k] -= BLOCK_SIZE * BLOCK_SIZE * d_f[k];
        tmp = d_newpressure;
        d_newpressure = d_oldpressure;
        d_oldpressure  = tmp;
    }
    d_pressure[k] = d_oldpressure[k];
    //if(d_abs(k_z - GRID_COUNT/2) * d_abs(k_z - GRID_COUNT/2) + 
    //d_abs(k_y - GRID_COUNT/2) * d_abs(k_y - GRID_COUNT/2) +
    //d_abs(k_x - GRID_COUNT/2) * d_abs(k_x - GRID_COUNT/2) < GRID_COUNT  *GRID_COUNT / (5*5*25)){
    //    printf("%f\n", d_pressure[k]);
    //}
}
__global__ void resetPressure(float* d_pressure){
    const int k_x = threadIdx.x + blockDim.x * blockIdx.x;
    const int k_y = threadIdx.y + blockDim.y * blockIdx.y;
    const int k_z = threadIdx.z + blockDim.z * blockIdx.z;
    if ((k_x >= dev_Ld[0] ) || (k_y >= dev_Ld[1] ) || (k_z >= dev_Ld[2])) return;
    const int k = flatten(k_x, k_y, k_z);
    d_pressure[k] = 0;
}
__global__ void substractPressureGradient(float3 * d_vel, float* d_pressure){
    const int k_x = threadIdx.x + blockDim.x * blockIdx.x;
    const int k_y = threadIdx.y + blockDim.y * blockIdx.y;
    const int k_z = threadIdx.z + blockDim.z * blockIdx.z;
    if ((k_x >= dev_Ld[0] ) || (k_y >= dev_Ld[1] ) || (k_z >= dev_Ld[2])) return;
    if(k_x == 0 || k_x >= GRID_COUNT-1  || k_y == 0 || k_y >= GRID_COUNT-1 || k_z == 0 || k_z >= GRID_COUNT-1) return;
    const int k = flatten(k_x, k_y, k_z);
    d_vel[k].x -= dev_Deltat[0] *  (d_pressure[flatten(k_x+1,k_y,k_z)] - d_pressure[k]) / BLOCK_SIZE;
    d_vel[k].y -= dev_Deltat[0] *  (d_pressure[flatten(k_x,k_y+1,k_z)] - d_pressure[k]) / BLOCK_SIZE;
    d_vel[k].z -= dev_Deltat[0] *  (d_pressure[flatten(k_x,k_y,k_z+1)] - d_pressure[k]) / BLOCK_SIZE;
}
void forceIncompressibility(float3 * d_vel, float* d_pressure){
    // TODO: CHOLESKI PREPROCESS
    const int NFLAT = GRID_COUNT * GRID_COUNT * GRID_COUNT;
    const dim3 gridSize(blocksNeeded(GRID_COUNT, M_i.x), 
                        blocksNeeded(GRID_COUNT, M_i.y), 
                        blocksNeeded(GRID_COUNT, M_i.z));
    float * d_f;
    HANDLE_ERROR(cudaMalloc(&d_f, NFLAT*sizeof(float)));
    float * d_pressure1;    
    HANDLE_ERROR(cudaMalloc(&d_pressure1, NFLAT*sizeof(float)));
    
    resetPressure<<<gridSize, M_i>>>(d_pressure);
    HANDLE_ERROR(cudaPeekAtLastError());
    HANDLE_ERROR(cudaDeviceSynchronize());
    
    prepareJacobiMethod<<<gridSize, M_i>>>(d_vel, d_f);
    HANDLE_ERROR(cudaPeekAtLastError());
    HANDLE_ERROR(cudaDeviceSynchronize());
    jacobiIterations<<<gridSize, M_i>>>(d_pressure, d_pressure1, d_f);
    HANDLE_ERROR(cudaPeekAtLastError());    
    HANDLE_ERROR(cudaDeviceSynchronize());
    substractPressureGradient<<<gridSize, M_i>>>(d_vel, d_pressure);
    HANDLE_ERROR(cudaPeekAtLastError());
    HANDLE_ERROR(cudaDeviceSynchronize());

    HANDLE_ERROR(cudaFree(d_f));
    HANDLE_ERROR(cudaFree(d_pressure1));
}