// --------------------------------------------------------------------------
// gMini,
// a minimal Glut/OpenGL app to extend                              
//
// Copyright(C) 2007-2009                
// Tamy Boubekeur
//                                                                            
// All rights reserved.                                                       
//                                                                            
// This program is free software; you can redistribute it and/or modify       
// it under the terms of the GNU General Public License as published by       
// the Free Software Foundation; either version 2 of the License, or          
// (at your option) any later version.                                        
//                                                                            
// This program is distributed in the hope that it will be useful,            
// but WITHOUT ANY WARRANTY; without even the implied warranty of             
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              
// GNU General Public License (http://www.gnu.org/licenses/gpl.txt)           
// for more details.                                                          
//                                                                          
// --------------------------------------------------------------------------

#include "Camera.h"
#include <GL/glew.h>
#include <GL/glut.h>
#include <GL/glu.h>
#include <iostream>
#include <string>

// ---------------------------------------------
// BEGIN : Code from SGI
// ---------------------------------------------

/*
 * Pass the x and y coordinates of the last and current positions of
 * the mouse, scaled so they are from (-1.0 ... 1.0).
 *
 * The resulting rotation is returned as a quaternion rotation in the
 * first paramater.
 */
void
trackball(float q[4], float p1x, float p1y, float p2x, float p2y);

void
negate_quat(float *q, float *qn);

/*
 * Given two quaternions, add them together to get a third quaternion.
 * Adding quaternions to get a compound rotation is analagous to adding
 * translations to get a compound translation.  When incrementally
 * adding rotations, the first argument here should be the new
 * rotation, the second and third the total rotation (which will be
 * over-written with the resulting new total rotation).
 */
void
add_quats(float *q1, float *q2, float *dest);

/*
 * A useful function, builds a rotation matrix in Matrix based on
 * given quaternion.
 */
void
build_rotmatrix(float m[4][4], float q[4]);

/*
 * This function computes a quaternion based on an axis (defined by
 * the given vector) and an angle about which to rotate.  The angle is
 * expressed in radians.  The result is put into the third argument.
 */
void
axis_to_quat(float a[3], float phi, float q[4]);

// ---------------------------------------------
// END : Code from SGI
// ---------------------------------------------

using namespace std;

static int _spinning, _moving;
static int _beginu, _beginv;
static float _curquat[4];
static float _x, _y, _z;
static float __zoom;
static bool ini = false;

Camera::Camera () {
  fovAngle = 45.0;
  aspectRatio = 1.0;
  nearPlane = 0.1;
  farPlane = 10000.0;
  
  spinning = 0;
  moving = 0;
  beginu = 0;
  beginv = 0;
  
  trackball (curquat, 0.0, 0.0, 0.0, 0.0);
  x = y = z = 0.0;
  _zoom = 3.0;
  
  mouseRotatePressed = false;
  mouseMovePressed = false;
  mouseZoomPressed = false;
  lastX = 0;
  lastY = 0;
  lastZoom = 0;
}


void Camera::resize (int _W, int _H) {
  H = _H;
  W = _W;
  glViewport (0, 0, (GLint)W, (GLint)H);
  glMatrixMode (GL_PROJECTION);
  glLoadIdentity ();
  aspectRatio = static_cast<float>(W)/static_cast<float>(H);
  gluPerspective (fovAngle, aspectRatio, nearPlane, farPlane);
  glMatrixMode (GL_MODELVIEW);
}


void Camera::initPos () {
  if (!ini) {
  _spinning = spinning;
  _moving = moving;;
  _beginu = beginu;
  _beginv = beginv;
  _curquat[0] = curquat[0];
  _curquat[1] = curquat[1];
  _curquat[2] = curquat[2];
  _curquat[3] = curquat[3];
  _x = x;
  _y = y;
  _z = z;;
  __zoom = _zoom;
  ini = true;
  } else {
    spinning = _spinning;
    moving = _moving;;
    beginu = _beginu;
    beginv = _beginv;
    curquat[0] = _curquat[0];
    curquat[1] = _curquat[1];
    curquat[2] = _curquat[2];
    curquat[3] = _curquat[3];
    x = _x;
    y = _y;
    z = _z;;
    _zoom = __zoom;
  } 
}


void Camera::move (float dx, float dy, float dz) {
  x += dx;
  y += dy;
  z += dz;
}


void Camera::beginRotate (int u, int v) {
  beginu = u; 
  beginv = v;
  moving = 1;
  spinning = 0;
}


void Camera::rotate (int u, int v) {
  if (moving) {
    trackball(lastquat,
	      (2.0 * beginu - W) / W,
	      (H - 2.0 * beginv) / H,
	      (2.0 * u - W) / W,
	      (H - 2.0 * v) / H);
    beginu = u;
    beginv = v;
    spinning = 1;
    add_quats (lastquat, curquat, curquat);
  }
}


void Camera::endRotate () {
  moving = false;
}


void Camera::zoom (float z) {
  _zoom += z;
}


void Camera::apply () {
  glLoadIdentity();
  glTranslatef (x, y, z);
  GLfloat m[4][4]; 
  build_rotmatrix(m, curquat);
  glTranslatef (0.0, 0.0, -_zoom);
  glMultMatrixf(&m[0][0]);
}


void Camera::getPos (float & X, float & Y, float & Z) {
  GLfloat m[4][4]; 
  build_rotmatrix(m, curquat);
  float _x = -x;
  float _y = -y;
  float _z = -z + _zoom;
  X = m[0][0] * _x +  m[0][1] * _y +  m[0][2] * _z;
  Y = m[1][0] * _x +  m[1][1] * _y +  m[1][2] * _z;
  Z = m[2][0] * _x +  m[2][1] * _y +  m[2][2] * _z;
}

void Camera::handleMouseClickEvent (int button, int state, int x, int y) {
	if (state == GLUT_UP) {
        mouseMovePressed = false;
        mouseRotatePressed = false;
        mouseZoomPressed = false;
    } else {
        if (button == GLUT_LEFT_BUTTON) {
            beginRotate (x, y);
            mouseMovePressed = false;
            mouseRotatePressed = true;
            mouseZoomPressed = false;
        } else if (button == GLUT_RIGHT_BUTTON) {
            lastX = x;
            lastY = y;
            mouseMovePressed = true;
            mouseRotatePressed = false;
            mouseZoomPressed = false;
        } else if (button == GLUT_MIDDLE_BUTTON) {
            if (mouseZoomPressed == false) {
                lastZoom = y;
                mouseMovePressed = false;
                mouseRotatePressed = false;
                mouseZoomPressed = true;
            }
        }
    }
}

void Camera::handleMouseMoveEvent (int x, int y) {
	if (mouseRotatePressed == true) 
        rotate (x, y);
    else if (mouseMovePressed == true) {
        move ((x-lastX)/static_cast<float>(getScreenWidth ()), (lastY-y)/static_cast<float>(getScreenHeight ()), 0.0);
        lastX = x;
        lastY = y;
    }
    else if (mouseZoomPressed == true) {
        zoom (float (y-lastZoom)/getScreenHeight ());
        lastZoom = y;
    }
}

// ---------------------------------------------
// BEGIN : Code from SGI
// ---------------------------------------------

#include <cstdio>
/*
 * (c) Copyright 1993, 1994, Silicon Graphics, Inc.
 * ALL RIGHTS RESERVED
 * Permission to use, copy, modify, and distribute this software for
 * any purpose and without fee is hereby granted, provided that the above
 * copyright notice appear in all copies and that both the copyright notice
 * and this permission notice appear in supporting documentation, and that
 * the name of Silicon Graphics, Inc. not be used in advertising
 * or publicity pertaining to distribution of the software without specific,
 * written prior permission.
 *
 * THE MATERIAL EMBODIED ON THIS SOFTWARE IS PROVIDED TO YOU "AS-IS"
 * AND WITHOUT WARRANTY OF ANY KIND, EXPRESS, IMPLIED OR OTHERWISE,
 * INCLUDING WITHOUT LIMITATION, ANY WARRANTY OF MERCHANTABILITY OR
 * FITNESS FOR A PARTICULAR PURPOSE.  IN NO EVENT SHALL SILICON
 * GRAPHICS, INC.  BE LIABLE TO YOU OR ANYONE ELSE FOR ANY DIRECT,
 * SPECIAL, INCIDENTAL, INDIRECT OR CONSEQUENTIAL DAMAGES OF ANY
 * KIND, OR ANY DAMAGES WHATSOEVER, INCLUDING WITHOUT LIMITATION,
 * LOSS OF PROFIT, LOSS OF USE, SAVINGS OR REVENUE, OR THE CLAIMS OF
 * THIRD PARTIES, WHETHER OR NOT SILICON GRAPHICS, INC.  HAS BEEN
 * ADVISED OF THE POSSIBILITY OF SUCH LOSS, HOWEVER CAUSED AND ON
 * ANY THEORY OF LIABILITY, ARISING OUT OF OR IN CONNECTION WITH THE
 * POSSESSION, USE OR PERFORMANCE OF THIS SOFTWARE.
 *
 * US Government Users Restricted Rights
 * Use, duplication, or disclosure by the Government is subject to
 * restrictions set forth in FAR 52.227.19(c)(2) or subparagraph
 * (c)(1)(ii) of the Rights in Technical Data and Computer Software
 * clause at DFARS 252.227-7013 and/or in similar or successor
 * clauses in the FAR or the DOD or NASA FAR Supplement.
 * Unpublished-- rights reserved under the copyright laws of the
 * United States.  Contractor/manufacturer is Silicon Graphics,
 * Inc., 2011 N.  Shoreline Blvd., Mountain View, CA 94039-7311.
 *
 * OpenGL(TM) is a trademark of Silicon Graphics, Inc.
 */
/*
 * Trackball code:
 *
 * Implementation of a virtual trackball.
 * Implemented by Gavin Bell, lots of ideas from Thant Tessman and
 *   the August '88 issue of Siggraph's "Computer Graphics," pp. 121-129.
 *
 * Vector manip code:
 *
 * Original code from:
 * David M. Ciemiewicz, Mark Grossman, Henry Moreton, and Paul Haeberli
 *
 * Much mucking with by:
 * Gavin Bell
 */
#if defined(_WIN32)
#pragma warning (disable:4244)          /* disable bogus conversion warnings */
#endif
#include <cmath>

/*
 * This size should really be based on the distance from the center of
 * rotation to the point on the object underneath the mouse.  That
 * point would then track the mouse as closely as possible.  This is a
 * simple example, though, so that is left as an Exercise for the
 * Programmer.
 */
#define TRACKBALLSIZE  (0.8f)

/*
 * Local function prototypes (not defined in trackball.h)
 */
static float tb_project_to_sphere(float, float, float);
static void normalize_quat(float [4]);

void
vzero(float *v)
{
    v[0] = 0.0;
    v[1] = 0.0;
    v[2] = 0.0;
}

void
vset(float *v, float x, float y, float z)
{
    v[0] = x;
    v[1] = y;
    v[2] = z;
}

void
vsub(const float *src1, const float *src2, float *dst)
{
    dst[0] = src1[0] - src2[0];
    dst[1] = src1[1] - src2[1];
    dst[2] = src1[2] - src2[2];
}

void
vcopy(const float *v1, float *v2)
{
    register int i;
    for (i = 0 ; i < 3 ; i++)
        v2[i] = v1[i];
}

void
vcross(const float *v1, const float *v2, float *cross)
{
    float temp[3];

    temp[0] = (v1[1] * v2[2]) - (v1[2] * v2[1]);
    temp[1] = (v1[2] * v2[0]) - (v1[0] * v2[2]);
    temp[2] = (v1[0] * v2[1]) - (v1[1] * v2[0]);
    vcopy(temp, cross);
}

float
vlength(const float *v)
{
    return sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

void
vscale(float *v, float div)
{
    v[0] *= div;
    v[1] *= div;
    v[2] *= div;
}

void
vnormal(float *v)
{
    vscale(v,1.0/vlength(v));
}

float
vdot(const float *v1, const float *v2)
{
    return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
}

void
vadd(const float *src1, const float *src2, float *dst)
{
    dst[0] = src1[0] + src2[0];
    dst[1] = src1[1] + src2[1];
    dst[2] = src1[2] + src2[2];
}

/*
 * Ok, simulate a track-ball.  Project the points onto the virtual
 * trackball, then figure out the axis of rotation, which is the cross
 * product of P1 P2 and O P1 (O is the center of the ball, 0,0,0)
 * Note:  This is a deformed trackball-- is a trackball in the center,
 * but is deformed into a hyperbolic sheet of rotation away from the
 * center.  This particular function was chosen after trying out
 * several variations.
 *
 * It is assumed that the arguments to this routine are in the range
 * (-1.0 ... 1.0)
 */
void
trackball(float q[4], float p1x, float p1y, float p2x, float p2y)
{
    float a[3]; /* Axis of rotation */
    float phi;  /* how much to rotate about axis */
    float p1[3], p2[3], d[3];
    float t;

    if (p1x == p2x && p1y == p2y) {
        /* Zero rotation */
        vzero(q);
        q[3] = 1.0;
        return;
    }

    /*
     * First, figure out z-coordinates for projection of P1 and P2 to
     * deformed sphere
     */
    vset(p1,p1x,p1y,tb_project_to_sphere(TRACKBALLSIZE,p1x,p1y));
    vset(p2,p2x,p2y,tb_project_to_sphere(TRACKBALLSIZE,p2x,p2y));

    /*
     *  Now, we want the cross product of P1 and P2
     */
    vcross(p2,p1,a);

    /*
     *  Figure out how much to rotate around that axis.
     */
    vsub(p1,p2,d);
    t = vlength(d) / (2.0*TRACKBALLSIZE);

    /*
     * Avoid problems with out-of-control values...
     */
    if (t > 1.0) t = 1.0;
    if (t < -1.0) t = -1.0;
    phi = 2.0 * asin(t);

    axis_to_quat(a,phi,q);
}

/*
 *  Given an axis and angle, compute quaternion.
 */
void
axis_to_quat(float a[3], float phi, float q[4])
{
    vnormal(a);
    vcopy(a,q);
    vscale(q,sin(phi/2.0));
    q[3] = cos(phi/2.0);
}

/*
 * Project an x,y pair onto a sphere of radius r OR a hyperbolic sheet
 * if we are away from the center of the sphere.
 */
static float
tb_project_to_sphere(float r, float x, float y)
{
    float d, t, z;

    d = sqrt(x*x + y*y);
    if (d < r * 0.70710678118654752440) {    /* Inside sphere */
        z = sqrt(r*r - d*d);
    } else {           /* On hyperbola */
        t = r / 1.41421356237309504880;
        z = t*t / d;
    }
    return z;
}

/*
 * Given two rotations, e1 and e2, expressed as quaternion rotations,
 * figure out the equivalent single rotation and stuff it into dest.
 *
 * This routine also normalizes the result every RENORMCOUNT times it is
 * called, to keep error from creeping in.
 *
 * NOTE: This routine is written so that q1 or q2 may be the same
 * as dest (or each other).
 */

#define RENORMCOUNT 97

void
negate_quat(float q[4], float nq[4])
{
    nq[0] = -q[0];
    nq[1] = -q[1];
    nq[2] = -q[2];
    nq[3] = q[3];
}

void
add_quats(float q1[4], float q2[4], float dest[4])
{
    static int count=0;
    float t1[4], t2[4], t3[4];
    float tf[4];

#if 0
printf("q1 = %f %f %f %f\n", q1[0], q1[1], q1[2], q1[3]);
printf("q2 = %f %f %f %f\n", q2[0], q2[1], q2[2], q2[3]);
#endif

    vcopy(q1,t1);
    vscale(t1,q2[3]);

    vcopy(q2,t2);
    vscale(t2,q1[3]);

    vcross(q2,q1,t3);
    vadd(t1,t2,tf);
    vadd(t3,tf,tf);
    tf[3] = q1[3] * q2[3] - vdot(q1,q2);

#if 0
printf("tf = %f %f %f %f\n", tf[0], tf[1], tf[2], tf[3]);
#endif

    dest[0] = tf[0];
    dest[1] = tf[1];
    dest[2] = tf[2];
    dest[3] = tf[3];

    if (++count > RENORMCOUNT) {
        count = 0;
        normalize_quat(dest);
    }
}

/*
 * Quaternions always obey:  a^2 + b^2 + c^2 + d^2 = 1.0
 * If they don't add up to 1.0, dividing by their magnitued will
 * renormalize them.
 *
 * Note: See the following for more information on quaternions:
 *
 * - Shoemake, K., Animating rotation with quaternion curves, Computer
 *   Graphics 19, No 3 (Proc. SIGGRAPH'85), 245-254, 1985.
 * - Pletinckx, D., Quaternion calculus as a basic tool in computer
 *   graphics, The Visual Computer 5, 2-13, 1989.
 */
static void
normalize_quat(float q[4])
{
    int i;
    float mag;

    mag = sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3]);
    for (i = 0; i < 4; i++) q[i] /= mag;
}

/*
 * Build a rotation matrix, given a quaternion rotation.
 *
 */
void
build_rotmatrix(float m[4][4], float q[4])
{
    m[0][0] = 1.0 - 2.0 * (q[1] * q[1] + q[2] * q[2]);
    m[0][1] = 2.0 * (q[0] * q[1] - q[2] * q[3]);
    m[0][2] = 2.0 * (q[2] * q[0] + q[1] * q[3]);
    m[0][3] = 0.0;

    m[1][0] = 2.0 * (q[0] * q[1] + q[2] * q[3]);
    m[1][1]= 1.0 - 2.0 * (q[2] * q[2] + q[0] * q[0]);
    m[1][2] = 2.0 * (q[1] * q[2] - q[0] * q[3]);
    m[1][3] = 0.0;

    m[2][0] = 2.0 * (q[2] * q[0] - q[1] * q[3]);
    m[2][1] = 2.0 * (q[1] * q[2] + q[0] * q[3]);
    m[2][2] = 1.0 - 2.0 * (q[1] * q[1] + q[0] * q[0]);
    m[2][3] = 0.0;

    m[3][0] = 0.0;
    m[3][1] = 0.0;
    m[3][2] = 0.0;
    m[3][3] = 1.0;
}

// ---------------------------------------------
// BEGIN : Code from SGI
// ---------------------------------------------
