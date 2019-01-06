/* dev_R3grid.cu
 * R3 under discretization (discretize functor) to a grid
 * Ernest Yeung  ernestyalumni@gmail.com
 * 20160728
 */
#include "dev_R3grid.cuh"

__constant__ int dev_Ld[3];

__host__ dev_Grid3d::dev_Grid3d( dim3 Ld_in) : Ld(Ld_in)
{
    int nflat = this->NFLAT();
    HANDLE_ERROR(
        cudaMalloc((void**)&this->dev_temperature, nflat*sizeof(float)) 
    );
    HANDLE_ERROR(
        cudaMalloc((void**)&this->dev_velocity, nflat*sizeof(float3)) 
    );
    HANDLE_ERROR(
        cudaMalloc((void**)&this->dev_smokeDensity, nflat*sizeof(float)) 
    );
    HANDLE_ERROR(
        cudaMalloc((void**)&this->dev_smokeVoxelRadiance, nflat*sizeof(float))
    );
    HANDLE_ERROR(
        cudaMalloc((void**)&this->dev_smokeVoxelTransparency, nflat*sizeof(float))
    );
}

/*
__host__ dev_Grid3d::~dev_Grid3d() {
    HANDLE_ERROR( 
        cudaFree( this->dev_rho ) );
    HANDLE_ERROR(
        cudaFree( this->dev_E ) );
    HANDLE_ERROR(
        cudaFree( this->dev_u ) );
    
}
* */

__host__ int dev_Grid3d :: NFLAT() {
    return Ld.x*Ld.y*Ld.z;
}	



__device__ dev_block3d :: dev_block3d(unsigned int N_x, unsigned int N_y, unsigned int N_z ) :
    N_is {N_x,N_y,N_z} 
{}


__device__ int dev_block3d :: flatten(
                            int i_x, int i_y, int i_z) {
    return i_x+i_y*N_is[0]+i_z*N_is[0]*N_is[1];

}	
    
