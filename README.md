# IG3DA: Interactive Wood Combustion for Botanical Tree Models [2017 Pirk et al.] 

The compiled binary require CUDA compute capability 61

`DemoVideo.ogv` is a video example 

## Requirements
CUDA 9.0 

g++-5 (at least g++<6, otherwise some compilation problems may appear)

libcublas and libcusparse

## Compilation
First, you need to edit the Makefile.

set CC and CPP to g++-5

set NVCC path

set COMPUTE_CAPABILITY to the capability of your card [https://en.wikipedia.org/wiki/CUDA#GPUs_supported]

then `make -j`

## Options
You can change all the physics parameters in `./src/physics/physics.h`
## Usage
The red line shows the external force applied to the fluid at each frame

### Commands:
------------------
 w: Toggle wireframe mode

drag+left button: rotate model

drag+right button: move model

drag+middle button / drag + left button + right button : zoom
 
 g: Toggle debug grid display
 
 s: Toggle smoke and temperature sources
 
 p: Pause/Unpause the physics simulation
 
 r: Reset the state of the physics engine
 
 left,right,up,down: Add external force to the fluid in the ground plane
 
 PageUp, PageDown: Add external force to the fluid along the vertical axis
 
 q, esc: Quit