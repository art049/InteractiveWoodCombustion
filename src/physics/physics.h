#pragma once

#include <cuda.h>
#include <cuda_runtime.h>
#include <cuda_gl_interop.h>
 

#include <GL/glut.h>

#include "../cuda_common/errors.h"
#include "../cuda_common/tex_anim2d.cuh"

#include "constants.h"
#include "dev_R3grid.cuh"
#include "R3grid.h"
#include "heat_3d.cuh"

extern dev_Grid3d dev_grid3d;
extern GPUAnim2dTex* testGPUAnim2dTex;

class Physics
{
private:
    dim3 dev_L3;
    dev_Grid3d * dev_grid3d;
public:
    Physics();
    ~Physics();
    void update();
};