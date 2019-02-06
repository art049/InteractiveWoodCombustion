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

static const float Deltat[1] {0.03f}; 
static const uint GRID_COUNT =  100;
static const float GRID_SIZE = 1;
static const float BLOCK_SIZE = GRID_SIZE/GRID_COUNT;
static const dim3 M_i { 8 , 8 , 8  };
static const float T_AMBIANT = 20.0f;
static const float P_ATM = 0.0f;
static const float BUOY_ALPHA = 3.3; // SMOKE DENSITY
static const float BUOY_BETA = 2.; // TEMPERATURE
static const uint SEMILAGRANGIAN_ITERS = 15;
static const float VORTICITY_EPSILON = 1e-6;
static const float TEMPERATURE_ALPHA = 8e-5;
static const float TEMPERATURE_GAMMA = -8e-7;
static const int PRESSURE_JACOBI_ITERATIONS = 5;
static const float SMOKE_EXTINCTION_COEFF = 8e1;
static const float SMOKE_ALBEDO = 1;
static const int SMOKE_CIRCULAR_RAY_COUNT = 200 ;
static const float SMOKE_LIGHT_ANGLE = M_PI / 6;
static const float SMOKE_RAY_DELTA_ANGLE = SMOKE_LIGHT_ANGLE/SMOKE_CIRCULAR_RAY_COUNT;
static const float3 SMOKE_LIGHT_POS = {-1,0.5,2};
static const float3 SMOKE_LIGHT_DIR = {-1,0.5,2};

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
    float * smokeQuadsPositions;
    uint * smokeIndexes;
    float * smokeQuadsColors;
    GLuint smokeQuadVBO, smokeQuadIndexVBO;
    
    GLuint smokeColorBufferObj = 0;
    cudaGraphicsResource *cuda_smokeColorBufferObj_resource;

public:
    Physics();
    ~Physics();
    void initSmokeQuads();
    void renderSmokeQuads();
    void renderGrid();
    void renderLightRays();
    inline void toggleGrid() { gridEnabled = !gridEnabled; };
    void render();
    void update();
    void reset();
};

#endif

