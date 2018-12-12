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

#include "GLError.h"

void checkGLExceptions () {
    GLenum glerr = glGetError ();
    switch (glerr) {
        case GL_INVALID_ENUM:
        throw Exception ("Invalid Enum");
        break;
        case GL_INVALID_VALUE:
        throw Exception ("Invalid Value");
        break;
        case GL_INVALID_OPERATION:
        throw Exception ("Invalid Operation");
        break;

        case GL_STACK_OVERFLOW:
        throw Exception ("Stack Overflow");
        break;
        case GL_STACK_UNDERFLOW:
        throw Exception ("Stack Underflow");
        break;
        case GL_OUT_OF_MEMORY:
        throw Exception ("Out of Memory");
        break;
        case GL_TABLE_TOO_LARGE:
        throw Exception ("Table Too Large");
        break;
        case GL_NO_ERROR:
        break;
        default:
        throw Exception ("Other Error");
        break;
    }
}