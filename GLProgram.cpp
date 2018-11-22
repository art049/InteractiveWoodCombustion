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

#include "GLProgram.h"

#include "GLError.h"

using namespace std;

GLProgram::GLProgram (const std::string & name) :
_id (glCreateProgram ()),
_name (name) {}

GLProgram::~GLProgram () {
    glDeleteProgram (_id);
}

void GLProgram::attach (GLShader * shader) {
    glAttachShader (_id, shader->id());
    _shaders.push_back (shader);
}

void GLProgram::detach (GLShader * shader) {
    for (unsigned int i = 0; i < _shaders.size (); i++)
        if (_shaders[i]->id () == shader->id ())
            glDetachShader (_id, shader->id());
}

void GLProgram::link () {
    glLinkProgram (_id);
    printOpenGLError ("Linking Program " + name ());
    GLint linked;
    glGetProgramiv (_id, GL_LINK_STATUS, &linked);
    if (!linked)
        throw Exception ("Shaders not linked: " + infoLog ());
}

void GLProgram::use () {
    glUseProgram (_id);
}

void GLProgram::stop () {
    glUseProgram (0);
}

std::string GLProgram::infoLog () {
    std::string infoLogStr = "";
    int infologLength = 0;
    glGetProgramiv (_id, GL_INFO_LOG_LENGTH, &infologLength);
    printOpenGLError ("Gathering Shader InfoLog for Program " + name ());
    if (infologLength > 0) {
        GLchar *str = new GLchar[infologLength];
        int charsWritten  = 0;
        glGetProgramInfoLog (_id, infologLength, &charsWritten, str);
        printOpenGLError ("Gathering Shader InfoLog for Program " + name ());
        infoLogStr  = std::string (str);
        delete [] str;
    }
    return infoLogStr;
}

GLint GLProgram::getUniformLocation (const std::string & uniformName) {
    const GLchar * cname = uniformName.c_str ();
    GLint loc = glGetUniformLocation (_id, cname);
    if (loc == -1)
        throw Exception (std::string ("Program Error: No such uniform named ") + uniformName);
    printOpenGLError ("Wrong Uniform Variable [" + uniformName + "] for Program [" + name () + "]");
    return loc;
}

void GLProgram::setUniform1f (GLint location, float value) {
    use ();
    glUniform1f (location, value);
}

void GLProgram::setUniform1f (const std::string & name, float value) {
    use ();
    glUniform1f (getUniformLocation (name), value);
}

void GLProgram::setUniform2f (GLint location, float value0, float value1) {
    use ();
    glUniform2f (location, value0, value1);
}

void GLProgram::setUniform2f (const std::string & name, float value0, float value1) {
    use ();
    glUniform2f (getUniformLocation (name), value0, value1);
}

void GLProgram::setUniform3f (GLint location, float value0, float value1, float value2) {
    use ();
    glUniform3f (location, value0, value1, value2);
}

void GLProgram::setUniform3f (const std::string & name, float value0, float value1, float value2) {
    use ();
    glUniform3f (getUniformLocation (name), value0, value1, value2);
}

void GLProgram::setUniform4f (GLint location, float value0, float value1, float value2, float value3) {
    use ();
    glUniform4f (location, value0, value1, value2, value3);
}

void GLProgram::setUniform4f (const std::string & name, float value0, float value1, float value2, float value3) {
    use ();
    glUniform4f (getUniformLocation (name), value0, value1, value2, value3);
}

void GLProgram::setUniformMatrix4fv (GLint location, const float * values) {
    use ();
    glUniformMatrix4fv (location, 1, GL_FALSE, values);
}

void GLProgram::setUniformMatrix4fv (const std::string & name, const float * values) {
    use ();
    setUniformMatrix4fv (getUniformLocation (name), values);
}

void GLProgram::setUniformNf (GLint location, unsigned int numValues, const float * values) {
    use ();
    switch (numValues) {
        case 1: glUniform1f (location, values[0]); break;
        case 2: glUniform2f (location, values[0], values[1]); break;
        case 3: glUniform3f (location, values[0], values[1], values[2]); break;
        case 4: glUniform4f (location, values[0], values[1], values[2], values[3]); break;
        default: throw Exception ("Program Error: Wrong number of values to set for uniform float array."); break;
    }
}

void GLProgram::setUniformNf (const std::string & name, unsigned int numValues, const float * values) {
    use ();
    GLint loc = getUniformLocation (name);
    switch (numValues) {
        case 1: glUniform1f (loc, values[0]); break;
        case 2: glUniform2f (loc, values[0], values[1]); break;
        case 3: glUniform3f (loc, values[0], values[1], values[2]); break;
        case 4: glUniform4f (loc, values[0], values[1], values[2], values[3]); break;
        default: throw Exception ("Wrong number of values to set for uniform float array " + name + "."); break;
    }
}

void GLProgram::setUniform1i (GLint location, int value) {
    use ();
    glUniform1i (location, value);
}

void GLProgram::setUniform1i (const std::string & name, int value) {
    use ();
    glUniform1i (getUniformLocation (name), value);
}

void GLProgram::setUniformNi (GLint location, unsigned int numValues, const int * values) {
    use ();
    switch (numValues) {
        case 1: glUniform1i (location, values[0]); break;
        case 2: glUniform2i (location, values[0], values[1]); break;
        case 3: glUniform3i (location, values[0], values[1], values[2]); break;
        case 4: glUniform4i (location, values[0], values[1], values[2], values[3]); break;
        default: throw Exception ("Program Error: Wrong number of values to set for uniform int array."); break;
    }
}

void GLProgram::setUniformNi (const std::string & name, unsigned int numValues, const int * values) {
    use ();
    GLint loc = getUniformLocation (name);
    switch (numValues) {
        case 1: glUniform1i (loc, values[0]); break;
        case 2: glUniform2i (loc, values[0], values[1]); break;
        case 3: glUniform3i (loc, values[0], values[1], values[2]); break;
        case 4: glUniform4i (loc, values[0], values[1], values[2], values[3]); break;
        default: throw Exception ("Program Error: Wrong number of values to set for uniform int array " + name + "."); break;
    }
}

void GLProgram::reload () {
    for (unsigned int i = 0; i < _shaders.size (); i++) {
        _shaders[i]->reload ();
        attach (_shaders[i]);
    }
    link ();
}

GLProgram * GLProgram::genVFProgram (const std::string & name,
                                     const std::string & vertexShaderFilename,
                                     const std::string & fragmentShaderFilename) {
    GLProgram * p = new GLProgram (name);
    GLShader * vs = new GLShader (name + " Vertex Shader", GL_VERTEX_SHADER);
    GLShader * fs = new GLShader (name + " Fragment Shader",GL_FRAGMENT_SHADER);
    vs->loadFromFile (vertexShaderFilename);
    vs->compile ();
    p->attach(vs);
    fs->loadFromFile (fragmentShaderFilename);
    fs->compile ();
    p->attach(fs);
    p->link();
    return p;
}
