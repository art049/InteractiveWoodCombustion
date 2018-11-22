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

#include <cstdio>
#include <cstdlib>

#include "Exception.h"

#define printOpenGLError(X) printOglError ((X), __FILE__, __LINE__)

/// Returns 1 if an OpenGL error occurred, 0 otherwise.
inline int printOglError (const std::string & msg, const char * file, int line) {
    GLenum glErr;
    int    retCode = 0;
    glErr = glGetError ();
    while (glErr != GL_NO_ERROR) {
        printf ("glError in file %s @ line %d: %s - %s\n", file, line, gluErrorString(glErr), msg.c_str ());
        retCode = 1;
        glErr = glGetError ();
    }
    return retCode;
}

/// Throws an expection if an error occurred.
void checkGLExceptions ();