/* tex_anim2d.h
 * 2-dim. GPU texture animation 
 * Ernest Yeung  ernestyalumni@gmail.com
 * 20160720
 */
#ifndef __TEX_ANIM2D_H__
#define __TEX_ANIM2D_H__

#include "../physics/heat_3d.cuh"  // BC
#include <stdio.h>

#define GL_GLEXT_PROTOTYPES // needed for identifier glGenBuffer, glBindBuffer, glBufferData, glDeleteBuffers

#include <GL/glut.h>

#include <cuda_runtime.h>
#include <cuda_gl_interop.h>

#include "errors.h"

#include <cstdlib> // exit

#define MAX(x, y) (((x) > (y)) ? (x) : (y) )
#define W 480
#define H 480
#define DT 1.f // source intensity increment


extern BC bc ; // Boundary conds
extern unsigned int slice;
extern uint const GRID_DEPTH;


struct GPUAnim2dTex {
	GLuint pixbufferObj ; // OpenGL pixel buffer object
	GLuint texObj       ; // OpenGL texture object	

	cudaGraphicsResource *cuda_pixbufferObj_resource;
 
	int width, height;
 
	GPUAnim2dTex( int w, int h ) {
		width  = w;
		height = h;

		pixbufferObj = 0 ;
		texObj       = 0 ;
	}

	~GPUAnim2dTex() {
		exitfunc();
	}

	void initGLUT(int *argc, char **argv) {
		glutInit(argc, argv);
		glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);
		glutInitWindowSize(width, height);
		glutCreateWindow("Temp. Vis.");

		gluOrtho2D(0, width, height, 0);
	}
	
	void initPixelBuffer() {
		glGenBuffers(1, &pixbufferObj );
		glBindBuffer(GL_PIXEL_UNPACK_BUFFER, pixbufferObj);
		glBufferData(GL_PIXEL_UNPACK_BUFFER, width*height*sizeof(GLubyte)*4, 0, 
			GL_STREAM_DRAW);
		glGenTextures(1, &texObj);
		glBindTexture(GL_TEXTURE_2D, texObj );
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		cudaGraphicsGLRegisterBuffer( &cuda_pixbufferObj_resource, pixbufferObj, 
			cudaGraphicsMapFlagsWriteDiscard);
	}
	
	void exitfunc() {
		if (pixbufferObj) {
			cudaGraphicsUnregisterResource( cuda_pixbufferObj_resource);
			glDeleteBuffers(1, &pixbufferObj);
			glDeleteTextures(1, &texObj);
		}
	}
};	
	

#endif // # __TEX_ANIM2D_H__ 
