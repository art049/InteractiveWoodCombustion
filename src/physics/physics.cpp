#include "physics.h"

#define GL_GLEXT_PROTOTYPES
Physics::Physics(){
    dev_L3.x = static_cast<unsigned int>(GRID_WIDTH);
    dev_L3.y = static_cast<unsigned int>(GRID_HEIGHT);
    dev_L3.z = static_cast<unsigned int>(GRID_DEPTH);
    dev_grid3d = new dev_Grid3d (dev_L3);
    // physics
    constexpr std::array<int,3> LdS {GRID_WIDTH, GRID_HEIGHT, GRID_DEPTH };
    constexpr std::array<float,3> ldS {1.f, 1.f, 1.f };
    HANDLE_ERROR(
        cudaMemcpyToSymbol( dev_Deltat, Deltat, sizeof(float)*1,0) );

    const float heat_params[2] { 
                                 0.00500f,
                                 1.f } ; // \kappa 
                                        // heat capacity for constant volume, per volume 

    HANDLE_ERROR(
        cudaMemcpyToSymbol( dev_heat_params, heat_params, sizeof(float)*2,0,cudaMemcpyHostToDevice) );
    
    const int Ld_to_const[3] { LdS[0], LdS[1], LdS[2] } ;
    
    HANDLE_ERROR(
        cudaMemcpyToSymbol( dev_Ld, Ld_to_const, sizeof(int)*3,0,cudaMemcpyHostToDevice) );
    
    Grid3d grid3d( LdS, ldS);
    
    const float hds[3] { grid3d.hd[0], grid3d.hd[1], grid3d.hd[2] } ;
    
    // sanity check
    //std::cout << " hds : .x : " << hds[0] << " .y : " << hds[1] << " .z : " << hds[2] << std::endl;
    
//	set1DerivativeParameters(hds);
    set2DerivativeParameters(hds);
//	set3DerivativeParameters(hds);
//	set4DerivativeParameters(hds);
    
    resetTemperature( dev_grid3d->dev_temperature, dev_L3, bc, M_i);
    
    testGPUAnim2dTex->initPixelBuffer();

}
Physics::~Physics() {
    HANDLE_ERROR(
        cudaFree(dev_grid3d->dev_temperature) 
    );
    HANDLE_ERROR(
        cudaFree(dev_grid3d->dev_velocity) 
    );
    HANDLE_ERROR(
        cudaFree(dev_grid3d->dev_smokeDensity)
    );
}
void Physics::update() {
    
    uchar4 *d_out = 0;
    cudaGraphicsMapResources(1, &testGPUAnim2dTex->cuda_pixbufferObj_resource, 0);
    cudaGraphicsResourceGetMappedPointer((void **)&d_out, NULL,
        testGPUAnim2dTex->cuda_pixbufferObj_resource);
    
    kernelLauncher(d_out, dev_grid3d->dev_temperature, dev_L3, bc, M_i, slice );

    cudaGraphicsUnmapResources(1, &testGPUAnim2dTex->cuda_pixbufferObj_resource, 0);
    
}