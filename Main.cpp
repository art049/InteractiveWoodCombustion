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

#include <GL/glew.h>
#include <GL/glut.h>
#include <iostream>
#include <vector>
#include <string>
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <cmath>

#include "Vec3.h"
#include "Camera.h"
#include "Mesh.h"
#include "GLProgram.h"
#include "Exception.h"

#include "LightSource.h"

using namespace std;

static const unsigned int DEFAULT_SCREENWIDTH = 1024;
static const unsigned int DEFAULT_SCREENHEIGHT = 768;
static const string DEFAULT_MESH_FILE ("models/sphere.off");

static const string appTitle ("IG3DA - Interactive Wood Combustion for Botanical Tree Models");
static const string myName ("Arthur Pastel");
static GLint window;
static unsigned int FPS = 0;
static bool fullScreen = false;

static Camera camera;
static Mesh mesh;
GLProgram * glProgram;

static std::vector<Vec3f> colorResponses; // Cached per-vertex color response, updated at each frame
static std::vector<LightSource> lightSources;
void printUsage () {
	std::cerr << std::endl 
                  << appTitle << std::endl
                  << "Author: " << myName << std::endl << std::endl
                  << "Usage: ./main [<file.off>]" << std::endl
                  << "Commands:" << std::endl 
                  << "------------------" << std::endl
                  << " ?: Print help" << std::endl
                  << " w: Toggle wireframe mode" << std::endl
                  << " <drag>+<left button>: rotate model" << std::endl 
                  << " <drag>+<right button>: move model" << std::endl
                  << " <drag>+<middle button>: zoom" << std::endl
                  << " q, <esc>: Quit" << std::endl << std::endl; 
}

void init (const char * modelFilename) {
    glewExperimental = GL_TRUE;
    glewInit (); // init glew, which takes in charges the modern OpenGL calls (v>1.2, shaders, etc)
    //glCullFace (GL_BACK);     // Specifies the faces to cull (here the ones pointing away from the camera)
    //glEnable (GL_CULL_FACE); // Enables face culling (based on the orientation defined by the CW/CCW enumeration).
    glDepthFunc (GL_LESS); // Specify the depth test for the z-buffer
    glEnable (GL_DEPTH_TEST); // Enable the z-buffer in the rasterization
    glEnableClientState (GL_VERTEX_ARRAY);
    glEnableClientState (GL_NORMAL_ARRAY);
    glEnableClientState (GL_COLOR_ARRAY);
    glEnable (GL_NORMALIZE);
	glLineWidth (2.0); // Set the width of edges in GL_LINE polygon mode
    glClearColor (0.0f, 0.0f, 0.0f, 1.0f); // Background color
	mesh.loadOFF (modelFilename);
    colorResponses.resize (mesh.positions ().size ());
    camera.resize (DEFAULT_SCREENWIDTH, DEFAULT_SCREENHEIGHT);
    try {
        glProgram = GLProgram::genVFProgram ("Simple GL Program", "shader.vert", "shader.frag"); // Load and compile pair of shaders
        glProgram->use (); // Activate the shader program

    } catch (Exception & e) {
        cerr << e.msg () << endl;
    }

    //Light config
    lightSources.push_back(LightSource(Vec3f(1,0,0),Vec3f(0,1,0),0.6));
    lightSources.push_back(LightSource(Vec3f(-1,0,0),Vec3f(1,0,0),0.6));
    lightSources.push_back(LightSource(Vec3f(0,0,-1),Vec3f(1,1,1),0.5));
    lightSources.push_back(LightSource(Vec3f(0,0,1),Vec3f(1,1,1),0.5));
}

// EXERCISE : the following color response shall be replaced with a proper reflectance evaluation/shadow test/etc.
void updatePerVertexColorResponse () {
  std::vector<Vec3f> positions = mesh.positions();
  std::vector<Vec3f> normals = mesh.normals();
  Vec3f campos;
  camera.getPos(campos);
  Vec3f lumpos;
  float intensity;
  float ks,kd;
  ks = 5.f;
  kd = 1.f;
  Vec3f wi,wh,wo,n,color,r;
  for (unsigned int i = 0; i < colorResponses.size (); i++){
    wo = campos - positions[i];
    wo.normalize();
    color = Vec3f(0,0,0);
    n = normals[i];
    for(vector<LightSource>::iterator it = lightSources.begin();it!=lightSources.end();it++){
      lumpos = it->getPos();
      intensity = it->getIntensityAt(positions[i]);
      wi = (positions[i]-lumpos);
      wi.normalize();
      wh = wi + wo;
      wh.normalize();
      //For phong
      r = 2.f*n*dot(wi,n) - wi;
      
      color +=(float) (intensity*(kd/M_PI+ks*pow(dot(n,wh),10.f)) * dot(n,wi))*it->getColor();
    }
    colorResponses[i] = color;
    
        
  }
}

void renderScene () {
    updatePerVertexColorResponse ();
    glVertexPointer (3, GL_FLOAT, sizeof (Vec3f), (GLvoid*)(&(mesh.positions()[0])));
    glNormalPointer (GL_FLOAT, 3*sizeof (float), (GLvoid*)&(mesh.normals()[0]));
    glColorPointer (3, GL_FLOAT, sizeof (Vec3f), (GLvoid*)(&(colorResponses[0])));
    glDrawElements (GL_TRIANGLES, 3*mesh.triangles().size(), GL_UNSIGNED_INT, (GLvoid*)((&mesh.triangles()[0])));
}

void reshape(int w, int h) {
    camera.resize (w, h);
}

void display () {
    glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    camera.apply (); 
    renderScene ();
    glFlush ();
    glutSwapBuffers (); 
}

void key (unsigned char keyPressed, int x, int y) {
    switch (keyPressed) {
    case 'f':
        if (fullScreen) {
            glutReshapeWindow (camera.getScreenWidth (), camera.getScreenHeight ());
            fullScreen = false;
        } else {
            glutFullScreen ();
            fullScreen = true;
        }      
        break;
    case 'q':
    case 27:
        exit (0);
        break;
    case 'w':
        GLint mode[2];
		glGetIntegerv (GL_POLYGON_MODE, mode);
		glPolygonMode (GL_FRONT_AND_BACK, mode[1] ==  GL_FILL ? GL_LINE : GL_FILL);
        break;
        break;
    default:
        printUsage ();
        break;
    }
}

void mouse (int button, int state, int x, int y) {
    camera.handleMouseClickEvent (button, state, x, y);
}

void motion (int x, int y) {
    camera.handleMouseMoveEvent (x, y);
}

void idle () {
    static float lastTime = glutGet ((GLenum)GLUT_ELAPSED_TIME);
    static unsigned int counter = 0;
    counter++;
    float currentTime = glutGet ((GLenum)GLUT_ELAPSED_TIME);
    if (currentTime - lastTime >= 1000.0f) {
        FPS = counter;
        counter = 0;
        static char winTitle [128];
        unsigned int numOfTriangles = mesh.triangles ().size ();
        sprintf (winTitle, "Number Of Triangles: %d - FPS: %d", numOfTriangles, FPS);
        string title = appTitle + " - By " + myName  + " - " + winTitle;
        glutSetWindowTitle (title.c_str ());
        lastTime = currentTime;
    }
    glutPostRedisplay (); 
}

int main (int argc, char ** argv) {
    if (argc > 2) {
        printUsage ();
        exit (1);
    }
    glutInit (&argc, argv);
    glutInitDisplayMode (GLUT_RGBA | GLUT_DEPTH | GLUT_DOUBLE);
    glutInitWindowSize (DEFAULT_SCREENWIDTH, DEFAULT_SCREENHEIGHT);
    window = glutCreateWindow (appTitle.c_str ());
    init (argc == 2 ? argv[1] : DEFAULT_MESH_FILE.c_str ());
    glutIdleFunc (idle);
    glutReshapeFunc (reshape);
    glutDisplayFunc (display);
    glutKeyboardFunc (key);
    glutMotionFunc (motion);
    glutMouseFunc (mouse);
    printUsage ();  
    glutMainLoop ();
    return 0;
}

