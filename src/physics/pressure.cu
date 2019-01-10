#include "pressure.cuh"

__global__ void prepareSystem(int NFLAT, float3* d_vel, float* d_b, float* d_val, int* d_cind, int * d_rptr) {
    const int k_x = threadIdx.x + blockDim.x * blockIdx.x;
    const int k_y = threadIdx.y + blockDim.y * blockIdx.y;
    const int k_z = threadIdx.z + blockDim.z * blockIdx.z;
    if ((k_x >= dev_Ld[0] ) || (k_y >= dev_Ld[1] ) || (k_z >= dev_Ld[2])) return;
    const int k = flatten(k_x, k_y, k_z);
    //B term
    d_b[k] = d_vel[flatten(k_x+1, k_y, k_z)].x - d_vel[flatten(k_x, k_y, k_z).x] + 
             d_vel[flatten(k_x, k_y+1, k_z)].y - d_vel[flatten(k_x, k_y, k_z).y] + 
             d_vel[flatten(k_x, k_y, k_z+1)].z - d_vel[flatten(k_x, k_y, k_z).z] ;
    d_b[k] /= BLOCK_SIZE * dev_Deltat[0];
    
    // Matrix 
    if(k_x > 0 && k_x < GRID_COUNT - 1 && k_y > 0 && k_y < GRID_COUNT - 1 && k_z > 0 && k_z < GRID_COUNT - 1) {
        const int offset = 7 * flatten(k_x-1, k_y-1, k_z-1, GRID_COUNT-2, GRID_COUNT-2, GRID_COUNT-2);
        d_val [offset    ] = -6; 
        d_cind[offset    ] = k;
        d_val [offset + 1] =  1; 
        d_cind[offset + 1] = flatten(k_x+1, k_y, k_z);
        d_val [offset + 2] =  1; 
        d_cind[offset + 1] = flatten(k_x-1, k_y, k_z);
        d_val [offset + 3] =  1; 
        d_cind[offset + 1] = flatten(k_x, k_y+1, k_z);
        d_val [offset + 4] =  1; 
        d_cind[offset + 1] = flatten(k_x, k_y-1, k_z);
        d_val [offset + 5] =  1; 
        d_cind[offset + 1] = flatten(k_x, k_y, k_z+1);
        d_val [offset + 6] =  1;
        d_cind[offset + 1] = flatten(k_x+1, k_y, k_z-1);
        
        d_rptr[k] 
    }
    //PRESSURE DIRICHLET BOUNDARY CONDITION
    else {

    }


}

void forceIncompressibility(float3 * d_vel, float* d_pressure){
    // TODO: CHOLESKI PREPROCESS
    const int NFLAT = GRID_COUNT * GRID_COUNT * GRID_COUNT;
    const dim3 gridSize(blocksNeeded(GRID_COUNT, M_i.x), 
                        blocksNeeded(GRID_COUNT, M_i.y), 
                        blocksNeeded(GRID_COUNT, M_i.z));

    // CGLS solver config
    float shift = 0;
    float tol = 1e-6;
    int maxit = 20;
    bool quiet = false;
    int m = NFLAT;
    int n = NFLAT;
    int nnz = 7 * NFLAT;
    float *d_val, *d_b;
    int *d_cind, *d_rptr;
    HANDLE_ERROR(cudaMalloc(&d_val, (nnz + m) * sizeof(float)));
    d_b = d_val + nnz;
    HANDLE_ERROR(cudaMalloc(&d_cind, (nnz + m + 1) * sizeof(int)));
    d_rptr = d_cind + nnz;

    prepareSystem<<<gridSize, M_i>>>(NFLAT, d_vel, d_b, d_val, d_cind, d_rptr);
    HANDLE_ERROR(cudaPeekAtLastError());

    int flag = cgls::Solve<float, cgls::CSR>(d_val, d_rptr, d_cind, m, n, nnz, d_b, d_pressure, 
                                             shift, tol, maxit, quiet);
    if (flag != 0)
        printf("[CGLS warning] Flag = %d, Error = %e\n", flag, err);
    HANDLE_ERROR(cudaFree(d_val));
    HANDLE_ERROR(cudaFree(d_cind));

    HANDLE_ERROR(cudaDeviceSynchronize());
}