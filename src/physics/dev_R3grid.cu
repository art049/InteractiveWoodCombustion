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
    int nfflat = this->NFFLAT();
    HANDLE_ERROR(
        cudaMalloc((void**)&this->dev_temperature0, nflat*sizeof(float)) 
    );
    HANDLE_ERROR(
        cudaMalloc((void**)&this->dev_temperature1, nflat*sizeof(float)) 
    );
    HANDLE_ERROR(
        cudaMalloc((void**)&this->dev_pressure, nflat*sizeof(float)) 
    );
    HANDLE_ERROR(
        cudaMalloc((void**)&this->dev_velocity0, nfflat*sizeof(float3)) 
    );
    HANDLE_ERROR(
        cudaMalloc((void**)&this->dev_velocity1, nfflat*sizeof(float3)) 
    );
    HANDLE_ERROR(
        cudaMalloc((void**)&this->dev_ccvelocity, nflat*sizeof(float3)) 
    );
    HANDLE_ERROR(
        cudaMalloc((void**)&this->dev_vorticity, nfflat*sizeof(float3)) //FFLAT ?
    );
    HANDLE_ERROR(
        cudaMalloc((void**)&this->dev_smokeDensity0, nflat*sizeof(float)) 
    );
    HANDLE_ERROR(
        cudaMalloc((void**)&this->dev_smokeDensity1, nflat*sizeof(float)) 
    );
    HANDLE_ERROR(
        cudaMalloc((void**)&this->dev_smokeVoxelRadiance, nflat*sizeof(float))
    );
    HANDLE_ERROR(
        cudaMalloc((void**)&this->dev_smokeVoxelTransparency, nflat*sizeof(float))
    );
}


__host__ dev_Grid3d::~dev_Grid3d() {
    HANDLE_ERROR(
        cudaFree(dev_temperature0) 
    );
    HANDLE_ERROR(
        cudaFree(dev_temperature1) 
    );
    HANDLE_ERROR(
        cudaFree(dev_pressure) 
    );
    HANDLE_ERROR(
        cudaFree(dev_velocity0) 
    );
    HANDLE_ERROR(
        cudaFree(dev_velocity1) 
    );
    HANDLE_ERROR(
        cudaFree(dev_ccvelocity) 
    );
    HANDLE_ERROR(
        cudaFree(dev_vorticity) 
    );
    HANDLE_ERROR(
        cudaFree(dev_smokeDensity0)
    );  
    HANDLE_ERROR(
        cudaFree(dev_smokeDensity1)
    );    
    HANDLE_ERROR(
        cudaFree(dev_smokeVoxelRadiance) 
    );
    HANDLE_ERROR(
        cudaFree(dev_smokeVoxelTransparency) 
    );
}


__host__ int dev_Grid3d :: NFLAT() {
    return Ld.x*Ld.y*Ld.z;
}	
__host__ int dev_Grid3d :: NFFLAT() {
    return (Ld.x+1)*(Ld.y+1)*(Ld.z+1);
}	



__device__ dev_block3d :: dev_block3d(unsigned int N_x, unsigned int N_y, unsigned int N_z ) :
    N_is {N_x,N_y,N_z} 
{}


__device__ int dev_block3d :: flatten(
                            int i_x, int i_y, int i_z) {
    return i_x+i_y*N_is[0]+i_z*N_is[0]*N_is[1];

}	
    
