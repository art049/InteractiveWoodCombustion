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

#include "GLShader.h"

#include <iostream>
#include <fstream>

#include "GLError.h"
#include "Exception.h"

using namespace std;

GLShader::GLShader (const std::string & name, GLuint type) {
	_id = glCreateShader (type);
	_name = name;
	_type = type;
	_filename = "";
	_source = "";
}

GLShader::~GLShader () {
	if (_id != 0)
		glDeleteShader (_id);
}

void GLShader::setSource (const std::string & source) {
	_source = source;
}

void GLShader::compile () {
	const GLchar * tmp = _source.c_str();
	glShaderSource (_id, 1, &tmp, NULL);
	glCompileShader (_id);
    printOpenGLError ("Compiling Shader " + name ());  // Check for OpenGL errors
    GLint shaderCompiled;
    glGetShaderiv (_id, GL_COMPILE_STATUS, &shaderCompiled);
    printOpenGLError ("Compiling Shader " + name ());  // Check for OpenGL errors
    if (!shaderCompiled)
      	throw Exception ("Error: shader not compiled. Info. Log.:\n" + infoLog () + "\nSource:\n" + _source);
}

void GLShader::loadFromFile (const std::string & filename) {
	_filename = filename;
	std::ifstream in (_filename.c_str ());
	if (!in)
		throw Exception ("Error loading shader source file: " + _filename);
	std::string source;
	char c[2];
	c[1]='\0';
	while (in.get (c[0]))
		source.append (c);
	in.close ();
	setSource (source);
}

void GLShader::reload () {
	if (_filename != "") {
		loadFromFile (std::string (_filename));
		compile ();
	}
}

std::string GLShader::infoLog () {
	std::string infoLogStr = "";
	int infologLength = 0;
	glGetShaderiv (_id, GL_INFO_LOG_LENGTH, &infologLength);
	printOpenGLError ("Gathering Shader InfoLog Length for " + name ());
	if (infologLength > 0) {
		GLchar *str = new GLchar[infologLength];
		int charsWritten  = 0;
		glGetShaderInfoLog (_id, infologLength, &charsWritten, str);
		printOpenGLError ("Gathering Shader InfoLog for " + name ());
		infoLogStr  = std::string (str);
		delete [] str;
	}
	return infoLogStr;
}