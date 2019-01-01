#ifndef __PHYSICS_H__
#define __PHYSICS_H__

#include <cuda.h>
#include <cuda_runtime.h>
#include <cuda_gl_interop.h>
 

#include <GL/glut.h>
#include "../cuda_common/errors.h"
#include "../cuda_common/tex_anim2d.cuh"


//#include "constants.h"
#include "dev_R3grid.cuh"
#include "R3grid.h"
#include "heat_3d.cuh"

//extern GPUAnim2dTex* testGPUAnim2dTex;

static const float Deltat[1] {0.0001f}; 
static const uint GRID_WIDTH =  100;
static const uint GRID_HEIGHT=  100;
static const uint GRID_DEPTH =  100;
static const dim3 M_i { 16 , 16 , 4  };


class Physics
{
private:
    dim3 dev_L3;
    dev_Grid3d * dev_grid3d;
    Grid3d * grid3d;
    float * smokeQuadsPositions;
    uint * smokeIndexes;
    float * smokeQuadsColors;

public:
    Physics();
    ~Physics();
    void initSmokeQuads();
    void renderSmokeQuads();
    void render();
    void update();
};

#endif