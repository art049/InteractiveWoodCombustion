/* tex_anim2d.cu
 * 2-dim. GPU texture animation 
 * Ernest Yeung  ernestyalumni@gmail.com
 * 20160720
 */
#include "tex_anim2d.cuh"
  
int iterationCount = 0 ;

// int x, y location of pipe center; float rad; radius of "pipe", int chamfer; float t_s, t_a, t_g
//BC bc = {W / 2, H / 2, W / 15.f, 150, 212.f, 70.f, 0.f}; // Boundary conds
BC bc = {W / 2, H / 2, W / 15.f, 150, 212.f, 10.f, 150.f}; // Boundary conds
unsigned int slice = GRID_COUNT/2;
// interactions

// make* functions make functions to pass into OpenGL (note OpenGL is inherently a C API

void make_draw_texture(int w, int h) {
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, w, h, 0, GL_RGBA, 
		GL_UNSIGNED_BYTE, NULL);
	glEnable(GL_TEXTURE_2D);
	glBegin(GL_QUADS);
	glTexCoord2f(0.0f, 0.0f); glVertex2f(0,0);
	glTexCoord2f(0.0f, 1.0f); glVertex2f(0,h);
	glTexCoord2f(1.0f, 1.0f); glVertex2f(w,h);
	glTexCoord2f(1.0f, 0.0f); glVertex2f(w,0);
	glEnd();
	glDisable(GL_TEXTURE_2D);
}	
