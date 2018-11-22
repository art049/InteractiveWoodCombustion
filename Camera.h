// --------------------------------------------------------------------------
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

#pragma once

#include "Vec3.h"

class Camera {
public:
  Camera ();
  virtual ~Camera () {}
  
  inline float getFovAngle () const { return fovAngle; }
  inline void setFovAngle (float newFovAngle) { fovAngle = newFovAngle; }
  inline float getAspectRatio () const { return aspectRatio; }
  inline float getNearPlane () const { return nearPlane; }
  inline void setNearPlane (float newNearPlane) { nearPlane = newNearPlane; }
  inline float getFarPlane () const { return farPlane; }
  inline void setFarPlane (float newFarPlane) { farPlane = newFarPlane; }
  inline unsigned int getScreenWidth () const { return W; }
  inline unsigned int getScreenHeight () const { return H; }
  
  void resize (int W, int H);
  
  void initPos ();

  void move (float dx, float dy, float dz);
  void beginRotate (int u, int v);
  void rotate (int u, int v);
  void endRotate ();
  void zoom (float z);
  void apply ();
  
  void getPos (float & x, float & y, float & z);
  inline void getPos (Vec3f & p) { getPos (p[0], p[1], p[2]); }
    
  // Connecting typical GLUT events
  void handleMouseClickEvent (int button, int state, int x, int y);
  void handleMouseMoveEvent (int x, int y);
  
private:
  float fovAngle;
  float aspectRatio;
  float nearPlane;
  float farPlane;
  
  int spinning, moving;
  int beginu, beginv;
  int H, W;
  float curquat[4];
  float lastquat[4];
  float x, y, z;
  float _zoom;
  
  bool mouseRotatePressed;
  bool mouseMovePressed;
  bool mouseZoomPressed;
  int lastX, lastY, lastZoom;
};

// Some Emacs-Hints -- please don't remove:
//
//  Local Variables:
//  mode:C++
//  tab-width:4
//  End:
