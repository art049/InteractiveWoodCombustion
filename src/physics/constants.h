#ifndef CONSTANTS_H
#define CONSTANTS_H
#include <cuda_runtime.h>

const float Deltat[1] {0.0001f}; // Easier to pass array to cuda
const uint GRID_WIDTH =  480;
const uint GRID_HEIGHT=  480;
const uint GRID_DEPTH =  288;
const dim3 M_i { 16 , 16 , 4  };

#endif