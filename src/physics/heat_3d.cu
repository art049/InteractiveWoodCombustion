/* heat_3d.cu
 * 3-dim. Laplace eq. (heat eq.) by finite difference with shared memory
 * Ernest Yeung  ernestyalumni@gmail.com
 * 20160729
 */
#include "heat_3d.cuh"

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
  return make_float3(a.x+b, a.y+b, a.z+b);
}
__device__ float3 operator*(const float &b, const float3 &a) {
    return make_float3(a.x+b, a.y+b, a.z+b);
  }
__global__ void resetKernel(float *d_temp, float3* d_vel, float* d_smokedensity, BC bc) {
    const int k_x = blockIdx.x*blockDim.x + threadIdx.x;
    const int k_y = blockIdx.y*blockDim.y + threadIdx.y;
    const int k_z = blockIdx.z*blockDim.z + threadIdx.z;

    if ((k_x >= dev_Ld[0]) || (k_y >= dev_Ld[1]) || (k_z >= dev_Ld[2])) return;
    d_temp[k_z*dev_Ld[0]*dev_Ld[1] + k_y*dev_Ld[0] + k_x] = T_AMBIANT;
    d_vel[k_z*dev_Ld[0]*dev_Ld[1] + k_y*dev_Ld[0] + k_x] = {0.f, 0.f, 0.f};
    d_smokedensity[k_z*dev_Ld[0]*dev_Ld[1] + k_y*dev_Ld[0] + k_x] = 0.f;
}


__global__ void tempKernel(float *d_temp, BC bc) {
    constexpr int NUS = 1;
    constexpr int radius = NUS;
    
    extern __shared__ float s_in[];
    // global indices
    const int k_x = threadIdx.x + blockDim.x * blockIdx.x;
    const int k_y = threadIdx.y + blockDim.y * blockIdx.y;
    const int k_z = threadIdx.z + blockDim.z * blockIdx.z;
    if ((k_x >= dev_Ld[0] ) || (k_y >= dev_Ld[1] ) || (k_z >= dev_Ld[2])) return;
    const int k = flatten(k_x, k_y, k_z, dev_Ld[0], dev_Ld[1],dev_Ld[2]);
    // local width and height
    const int3 S = { static_cast<int>(blockDim.x + 2 * radius), 
                        static_cast<int>(blockDim.y + 2 * radius), 
                        static_cast<int>(blockDim.z + 2 * radius) };

    // local indices
    const int s_x = threadIdx.x + radius;
    const int s_y = threadIdx.y + radius;
    const int s_z = threadIdx.z + radius;
    const int s_k = flatten(s_x, s_y, s_z, S.x, S.y, S.z);
    // assign default color values for d_out (black)

    // Load regular cells
    s_in[s_k] = d_temp[k];
    // Load halo cells
    if (threadIdx.x < radius ) {
        s_in[flatten(s_x - radius, s_y, s_z, S.x, S.y, S.z)] = 
            d_temp[flatten(k_x - radius, k_y, k_z, dev_Ld[0], dev_Ld[1], dev_Ld[2])];
        s_in[flatten(s_x + blockDim.x, s_y, s_z, S.x, S.y, S.z)] = 
            d_temp[flatten(k_x + blockDim.x, k_y, k_z, dev_Ld[0], dev_Ld[1], dev_Ld[2])];
    }
    if (threadIdx.y < radius) {
        s_in[flatten(s_x, s_y - radius, s_z,S.x, S.y, S.z)] = 
            d_temp[flatten(k_x, k_y - radius, k_z, dev_Ld[0], dev_Ld[1], dev_Ld[2])];
        s_in[flatten(s_x, s_y + blockDim.y, s_z, S.x, S.y, S.z)] = 
            d_temp[flatten(k_x, k_y + blockDim.y, k_z, dev_Ld[0], dev_Ld[1], dev_Ld[2])];
    }
    if (threadIdx.z < radius) {
        s_in[flatten(s_x, s_y, s_z - radius, S.x, S.y, S.z)] = 
            d_temp[flatten(k_x, k_y, k_z - radius, dev_Ld[0], dev_Ld[1], dev_Ld[2])];
        s_in[flatten(s_x, s_y, s_z + blockDim.z, S.x, S.y, S.z)] = 
            d_temp[flatten(k_x, k_y , k_z + blockDim.z, dev_Ld[0], dev_Ld[1], dev_Ld[2])];
    }

    
    // Calculate squared distance from pipe center
    float dSq = ((k_x - bc.x)*(k_x - bc.x) + (k_y - bc.y)*(k_y - bc.y) + 
                    (k_z - dev_Ld[2]/2)*(k_z - dev_Ld[2]/2)  // this can be changed manually, to place the "pipe" source in the "middle" of the z-axis
    );
    // If inside pipe, set temp to t_s and return
    if (dSq < bc.rad*bc.rad) {
        d_temp[k] = bc.t_s;
        return;
    }

    /*
    // If outside plate, set temp to t_a and return
    if ((k_x == 0 ) || (k_x == dev_Ld[0] - 1) || (k_y == 0 ) ||
        (k_x + k_y < bc.chamfer) || (k_x - k_y > dev_Ld[0] - bc.chamfer)) {
            d_temp[k] = bc.t_a;
            return;
    }*/
    // boundary conditions, BC, for "sides" 
    if ((k_y == 0) || (k_y == dev_Ld[1]-1) || (k_z == 0) || (k_z == dev_Ld[2]-1) ) {
        d_temp[k] = bc.t_a;
        return; 
    }

    /*
        // If point is below ground, set temp to t_g and return
        if (k_y == dev_Ld[1] - 1) {
            d_temp[k] = bc.t_g;
            return;
        }
    */
    // If point is in front of inlet, set temp to t_g and return
    if (k_x == 0) {
        d_temp[k] = bc.t_g;
        return;
    }

    __syncthreads();
    // For all the remaining points, find temperature.
    
    float3 stencil[NUS][2] ; 
    
    const float centerval { s_in[ s_k] };
    
    for (int nu = 0; nu < NUS; ++nu) {
        stencil[nu][0].x = s_in[flatten(s_x-(nu+1),s_y,s_z,S.x,S.y,S.z)] ; 
        stencil[nu][1].x = s_in[flatten(s_x+(nu+1),s_y,s_z,S.x,S.y,S.z)] ; 
        stencil[nu][0].y = s_in[flatten(s_x,s_y-(nu+1),s_z,S.x,S.y,S.z)] ; 
        stencil[nu][1].y = s_in[flatten(s_x,s_y+(nu+1),s_z,S.x,S.y,S.z)] ; 
        stencil[nu][0].z = s_in[flatten(s_x,s_y,s_z-(nu+1),S.x,S.y,S.z)] ; 
        stencil[nu][1].z = s_in[flatten(s_x,s_y,s_z+(nu+1),S.x,S.y,S.z)] ; 
    }

    float tempval { dev_lap1( centerval, stencil ) };

    __syncthreads();

    d_temp[k] += dev_Deltat[0]*(dev_heat_params[0]/dev_heat_params[1])*tempval;

}


__global__ void velocityKernel(float *d_temp, float3* d_vel, float* d_smokedensity){
    const int k_x = threadIdx.x + blockDim.x * blockIdx.x;
    const int k_y = threadIdx.y + blockDim.y * blockIdx.y;
    const int k_z = threadIdx.z + blockDim.z * blockIdx.z;
    if ((k_x >= dev_Ld[0] ) || (k_y >= dev_Ld[1] ) || (k_z >= dev_Ld[2])) return;
    const int k = flatten(k_x, k_y, k_z, dev_Ld[0], dev_Ld[1],dev_Ld[2]);
    // External forces
    float3 fext = {0,0,0};
    float3 fconf = {0,0,0};// TODO
    float3 fbuoy = {0,0, -BUOY_ALPHA*d_smokedensity[k] + BUOY_BETA*(d_temp[k] - T_AMBIANT)};
    float3 f = fext + fconf + fbuoy;
    d_vel[k] = d_vel[k] + f * dev_Deltat[0];//NEED AVERAGE TO CENTER Appendix A
    float3 pos = make_float3(k_x*BLOCK_SIZE, k_y*BLOCK_SIZE, k_z*BLOCK_SIZE);
    
    // Advection
    // Iteratively compute alpha_m
    float3 alpha_m = d_vel[k] * dev_Deltat[0];
    for(uint i = 0; i < SEMILAGRANGIAN_ITERS; i++){
        float3 estimated = pos - alpha_m;
        int3 b = {static_cast<int>(estimated.x/BLOCK_SIZE),
                  static_cast<int>(estimated.y/BLOCK_SIZE),
                  static_cast<int>(estimated.z/BLOCK_SIZE)};
        float3 localCoord = (estimated - make_float3(b.x*BLOCK_SIZE, b.y*BLOCK_SIZE, b.z*BLOCK_SIZE)) * (1 / BLOCK_SIZE);
        alpha_m = localCoord.x     * d_vel[flatten(b.x, b.y, b.z)  ]+
                  (1-localCoord.x) * d_vel[flatten(b.x+1, b.y, b.z)]+
                  localCoord.y     * d_vel[flatten(b.x, b.y, b.z)  ]+
                  (1-localCoord.y) * d_vel[flatten(b.x, b.y+1, b.z)]+
                  localCoord.z     * d_vel[flatten(b.x, b.y, b.z)  ]+
                  (1-localCoord.z) * d_vel[flatten(b.x, b.y, b.z+1)];
        alpha_m = alpha_m * dev_Deltat[0];
    }
    // Backtracing 
    float3 estimated = pos - 2 * alpha_m;
    int3 b = {static_cast<int>(estimated.x/BLOCK_SIZE),
              static_cast<int>(estimated.y/BLOCK_SIZE),
              static_cast<int>(estimated.z/BLOCK_SIZE)};
    float3 localCoord = (estimated - make_float3(b.x*BLOCK_SIZE, b.y*BLOCK_SIZE, b.z*BLOCK_SIZE)) * (1 / BLOCK_SIZE);
    float3 dv= localCoord.x     * d_vel[flatten(b.x, b.y, b.z)  ]+
               (1-localCoord.x) * d_vel[flatten(b.x+1, b.y, b.z)]+
               localCoord.y     * d_vel[flatten(b.x, b.y, b.z)  ]+
               (1-localCoord.y) * d_vel[flatten(b.x, b.y+1, b.z)]+
               localCoord.z     * d_vel[flatten(b.x, b.y, b.z)  ]+
               (1-localCoord.z) * d_vel[flatten(b.x, b.y, b.z+1)];
    dv = dv * 2 * dev_Deltat[0];
    //NEED OLD GRID FOR THIS
    d_vel[k] = d_vel[k] + dv;

}

__global__ void smokeAdvectionKernel(float *d_temp, float3* d_vel, float* d_smoke){
    const int k_x = threadIdx.x + blockDim.x * blockIdx.x;
    const int k_y = threadIdx.y + blockDim.y * blockIdx.y;
    const int k_z = threadIdx.z + blockDim.z * blockIdx.z;
    if ((k_x >= dev_Ld[0] ) || (k_y >= dev_Ld[1] ) || (k_z >= dev_Ld[2])) return;
    const int k = flatten(k_x, k_y, k_z, dev_Ld[0], dev_Ld[1],dev_Ld[2]);
    
    // Advection
    // Iteratively compute alpha_m
    float3 pos = make_float3(k_x*BLOCK_SIZE, k_y*BLOCK_SIZE, k_z*BLOCK_SIZE);
    float3 alpha_m = d_vel[k] * dev_Deltat[0];
    for(uint i = 0; i < SEMILAGRANGIAN_ITERS; i++){
        float3 estimated = pos - alpha_m;
        int3 b = {static_cast<int>(estimated.x/BLOCK_SIZE),
                  static_cast<int>(estimated.y/BLOCK_SIZE),
                  static_cast<int>(estimated.z/BLOCK_SIZE)};
        float3 localCoord = (estimated - make_float3(b.x*BLOCK_SIZE, b.y*BLOCK_SIZE, b.z*BLOCK_SIZE)) * (1 / BLOCK_SIZE);
        alpha_m = localCoord.x     * d_vel[flatten(b.x, b.y, b.z)  ]+
                  (1-localCoord.x) * d_vel[flatten(b.x+1, b.y, b.z)]+
                  localCoord.y     * d_vel[flatten(b.x, b.y, b.z)  ]+
                  (1-localCoord.y) * d_vel[flatten(b.x, b.y+1, b.z)]+
                  localCoord.z     * d_vel[flatten(b.x, b.y, b.z)  ]+
                  (1-localCoord.z) * d_vel[flatten(b.x, b.y, b.z+1)];
        alpha_m = alpha_m * dev_Deltat[0];
    }
    // Backtracing 
    float3 estimated = pos - 2 * alpha_m;
    int3 b = {static_cast<int>(estimated.x/BLOCK_SIZE),
              static_cast<int>(estimated.y/BLOCK_SIZE),
              static_cast<int>(estimated.z/BLOCK_SIZE)};
    float3 localCoord = (estimated - make_float3(b.x*BLOCK_SIZE, b.y*BLOCK_SIZE, b.z*BLOCK_SIZE)) * (1 / BLOCK_SIZE);
    float ds= localCoord.x     * d_smoke[flatten(b.x, b.y, b.z)  ]+
               (1-localCoord.x) * d_smoke[flatten(b.x+1, b.y, b.z)]+
               localCoord.y     * d_smoke[flatten(b.x, b.y, b.z)  ]+
               (1-localCoord.y) * d_smoke[flatten(b.x, b.y+1, b.z)]+
               localCoord.z     * d_smoke[flatten(b.x, b.y, b.z)  ]+
               (1-localCoord.z) * d_smoke[flatten(b.x, b.y, b.z+1)];
    ds = ds * 2 * dev_Deltat[0];
    //NEED OLD GRID FOR THIS
    d_smoke[k] = d_smoke[k] + ds;
}

__global__ void float_to_char( uchar4* dev_out, const float* outSrc, unsigned int slice = 5000) {
    const int k_x = threadIdx.x + blockDim.x * blockIdx.x;
    const int k_y = threadIdx.y + blockDim.y * blockIdx.y;
    
    // choose at which z coordinate to make the slice in x-y plane
    const int zcoordslice = slice == 5000 ? 140 : slice; 
    
    const int k   = k_x + k_y * blockDim.x*gridDim.x ; 
    const int fulloffset = k + zcoordslice*blockDim.x*gridDim.x*blockDim.y*gridDim.y; 

    dev_out[k].x = 0;
    dev_out[k].z = 0;
    dev_out[k].y = 0;
    dev_out[k].w = 255;


    const unsigned char intensity = clip((int) outSrc[fulloffset] ) ;
    dev_out[k].x = intensity ;       // higher temp -> more red
    dev_out[k].z = 255 - intensity ; // lower temp -> more blue
    
}
__global__ void generateSmokeColorBuffer( uchar4* dev_out, const float* d_smoke) {
    const int k_x = threadIdx.x + blockDim.x * blockIdx.x;
    const int k_y = threadIdx.y + blockDim.y * blockIdx.y;
    const int k_z = threadIdx.z + blockDim.z * blockIdx.z;
    if ((k_x >= dev_Ld[0] ) || (k_y >= dev_Ld[1] ) || (k_z >= dev_Ld[2])) return;
    const int k = flatten(k_x, k_y, k_z, dev_Ld[0], dev_Ld[1],dev_Ld[2]);
    const unsigned char intensity = clip((int) (d_smoke[k]*255));
    for(uint i = 0; i < 4; i++){
        dev_out[4*k+i].x = intensity;
        dev_out[4*k+i].z = intensity;
        dev_out[4*k+i].y = intensity;
        dev_out[4*k+i].w = 0;
    }
}


void kernelLauncher(uchar4 *d_out, float *d_temp, float3* d_vel, float* d_smokedensity, dim3 Ld, BC bc, dim3 M_in, unsigned int slice) {
    const dim3 gridSize(blocksNeeded(Ld.x, M_in.x), blocksNeeded(Ld.y, M_in.y), 
                        blocksNeeded(Ld.z,M_in.z));
    const size_t smSz = (M_in.x + 2 * RAD)*(M_in.y + 2 * RAD)*(M_in.z + 2 * RAD)*sizeof(float);

    velocityKernel<<<gridSize, M_in, smSz>>>(d_temp, d_vel, d_smokedensity);
    smokeAdvectionKernel<<<gridSize, M_in, smSz>>>(d_temp, d_vel, d_smokedensity);
    //generateSmokeColorBuffer<<<gridSize, M_in, smSz>>>(d_out, d_smokedensity);
    //tempKernel<<<gridSize, M_in, smSz>>>(d_temp, bc);
    
    const dim3 out_gridSize( gridSize.x, gridSize.y );
    const dim3 out_M( M_in.x, M_in.y );
//
    float_to_char<<<out_gridSize,out_M>>>(d_out, d_temp, slice) ; 
}

void resetVariables(float *d_temp, float3* d_vel, float* d_smokedensity, dim3 Ld, BC bc, dim3 M_in) {
    const dim3 gridSize( blocksNeeded(Ld.x, M_in.x), blocksNeeded( Ld.y, M_in.y), 
                            blocksNeeded(Ld.z, M_in.z));

    resetKernel<<<gridSize, M_in>>>(d_temp, d_vel, d_smokedensity, bc);
}