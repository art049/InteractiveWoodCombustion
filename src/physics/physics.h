#ifndef __PHYSICS_H__
#define __PHYSICS_H__
#include <cuda.h>
#include <cuda_runtime.h>
#include <cuda_gl_interop.h>
#include <math.h>
 

//#include <GL/glut.h>
#include "../cuda_common/errors.h"
#include "../cuda_common/tex_anim2d.cuh"


//#include "constants.h"
#include "dev_R3grid.cuh"
#include "R3grid.h"
#include "heat_3d.cuh"
#include "smoke_render.cuh"
//extern GPUAnim2dTex* testGPUAnim2dTex;

static const float Deltat[1] {0.010f}; 
static const uint GRID_COUNT =  10;
static const float GRID_SIZE = 1;
static const float BLOCK_SIZE = GRID_SIZE/GRID_COUNT;
static const dim3 M_i { 8 , 8 , 8  };
static const float T_AMBIANT = 20.0f;
static const float P_ATM = 0.0f;
static const float BUOY_ALPHA = 0.6; // SMOKE DENSITY
static const float BUOY_BETA = 0.1; // TEMPERATURE
static const uint SEMILAGRANGIAN_ITERS = 5;
static const float VORTICITY_EPSILON = 1e-1;
static const float TEMPERATURE_ALPHA = 8e-5;//8e-5;
static const float TEMPERATURE_GAMMA = -4e-7;//8e-7;
static const int PRESSURE_JACOBI_ITERATIONS = 20;
static const float SMOKE_EXTINCTION_COEFF = 15e1;
static const float SMOKE_ALBEDO = 0.4;
static const int SMOKE_RAY_SQRT_COUNT = 300 ;
static const float3 SMOKE_LIGHT_DIR = {1,0,0};
static const float3 SMOKE_LIGHT_POS = {-1,0,0};
static const float SMOKE_LIGHT_RADIANCE = 5e0;
static const float EXTERNAL_FORCE_DELTA = 0.1;
static const float heat_params[2] { 
                                 0.00500f,
                                 1.f } ; // \kappa 
                                        // heat capacity for constant volume, per volume 

class Physics
{
private:
    dim3 dev_L3;
    dev_Grid3d * dev_grid3d;
    Grid3d * grid3d;
    int activeBuffer = 0;
    bool gridEnabled = false;
    bool raysEnabled = false;
    bool sourcesEnabled = true;
    float * smokeQuadsPositions;
    uint * smokeIndexes;
    float * smokeQuadsColors;
    GLuint smokeQuadVBO;
    GLuint smokeQuadIndexBO;
    float3 externalForce;
    GLuint smokeColorBufferObj = 0;
    cudaGraphicsResource *cuda_smokeColorBufferObj_resource;

public:
    Physics();
    ~Physics();
    void initSmokeQuads();
    void renderSmokeQuads(uint cameraAxis);
    void renderGrid();
    void renderLightRays();
    void renderExternalForce();
    inline void toggleGrid() { gridEnabled = !gridEnabled; };
    inline void toggleSources() { sourcesEnabled = !sourcesEnabled; };
    void render(uint cameraAxis);
    void update();
    void reset();
    void addExternalForce(float3 f);
};

#endif

