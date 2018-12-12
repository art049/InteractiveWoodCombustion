// --------------------------------------------------------------------------
// Copyright(C) 2009-2016
// Tamy Boubekeur
// 
// Permission granted to use this code only for teaching projects and 
// private practice.
//
// Do not distribute this code outside the teaching assignements.                                                                           
// All rights reserved.                                                       
// --------------------------------------------------------------------------
#pragma once

#include <GL/glew.h>

#include <string>
#include <vector>

#include "Exception.h"
#include "GLShader.h"

class GLProgram {
public:
  GLProgram (const std::string & name);
  virtual ~GLProgram ();

  inline GLuint id () const { return _id; }
  std::string name () const { return _name; }
  void attach (GLShader * shader);
  void detach (GLShader * shader);
  void link ();
  void use ();
  static void stop ();
  GLint getUniformLocation (const std::string & uniformName);
  void setUniform1f (GLint location, float value);
  void setUniform1f (const std::string & name, float value);
  void setUniform2f (GLint location, float value0, float value1);
  void setUniform2f (const std::string & name, float value0, float value1);
  void setUniform3f (GLint location, float value0, float value1, float value2);
  void setUniform3f (const std::string & name, float value0, float value1, float vlaue2);
  void setUniform4f (GLint location, float value0, float value1, float value2, float value3);
  void setUniform4f (const std::string & name, float value0, float value1, float value2, float value3);
  void setUniformMatrix4fv (GLint location, const float * values);
  void setUniformMatrix4fv (const std::string & name, const float * values);
  void setUniformNf (GLint location, unsigned int numValues, const float * values);
  void setUniformNf (const std::string & name, unsigned int numValues, const float * values);
  void setUniform1i (GLint location, int value);
  void setUniform1i (const std::string & name, int value);
  void setUniformNi (GLint location, unsigned int numValues, const int * values);
  void setUniformNi (const std::string & name, unsigned int numValues, const int * values);
  void reload ();
  
  // generate a simple program, with only vertex and fragment shaders.
  static GLProgram * genVFProgram (const std::string & name,
				                         const std::string & vertexShaderFilename,
                        				 const std::string & fragmentShaderFilename);

protected:
  std::string infoLog ();

private:
  GLuint _id;
  std::string _name;
  std::vector<GLShader*>_shaders;
};
