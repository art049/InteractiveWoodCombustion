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

class GLShader {
public:
  GLShader (const std::string & name, GLuint _type);
  virtual ~GLShader ();
  inline GLuint id () const { return _id; }
  inline const std::string & name () const { return _name; }
  inline GLenum type () const { return _type; }
  inline const std::string & source () const { return _source; }
  inline const std::string & filename () const { return _filename; }
  void setSource (const std::string & source);
  void compile ();
  void loadFromFile (const std::string & filename);
  void reload ();

 protected:
  std::string infoLog ();

 private:
  GLuint _id;
  std::string _name;
  GLuint _type;
  std::string _filename;
  std::string _source;
};