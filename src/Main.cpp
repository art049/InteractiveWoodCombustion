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
#include "TreeGraph.h"
#include "Tree.h"
#include "GLProgram.h"
#include "Exception.h"

#include "LightSource.h"

using namespace std;

static const unsigned int DEFAULT_SCREENWIDTH = 2048;
static const unsigned int DEFAULT_SCREENHEIGHT = 1536;
static unsigned int WIDTH = DEFAULT_SCREENWIDTH;
static unsigned int HEIGHT = DEFAULT_SCREENHEIGHT;

static const string appTitle ("Interactive Wood Combustion for Botanical Tree Models");
static const string myName ("Arthur Pastel");
static GLint window;
static unsigned int FPS = 0;
static bool fullScreen = false;
static GLuint smokeTexture;
static unsigned char * smokeImage;
static Camera camera;
static TreeGraph tree_graph(FIVE_BRANCH);
static Tree tree(&tree_graph);
GLProgram * glProgram;
static Vec3f frustrumRays[4], initialCampos;
static std::vector<Vec3f> colorResponses; // Cached per-vertex color response, updated at each frame
static std::vector<LightSource> lightSources;
void printUsage () {
	std::cerr << std::endl 
                  << appTitle << std::endl
                  << "Author: " << myName << std::endl << std::endl
                  << "Commands:" << std::endl 
                  << "------------------" << std::endl
                  << " ?: Print help" << std::endl
                  << " w: Toggle wireframe mode" << std::endl
                  << " <drag>+<left button>: rotate model" << std::endl 
                  << " <drag>+<right button>: move model" << std::endl
                  << " <drag>+<middle button>: zoom" << std::endl
                  << " q, <esc>: Quit" << std::endl << std::endl; 
}
void genCheckerboard (unsigned int width, unsigned int height, unsigned char * image){
    const int psize = 20;
    const int fixedAlpha = 150;
    bool oddh,oddw;
    unsigned int index = 0;
    for(unsigned int y=0;y<height;y++) {
        oddh = (y / psize) % 2 == 1;
        for (unsigned int x = 0; x < width; x++) {
            oddw = (x / psize) % 2 == 1;
            if (oddh == oddw) {
                image[index++] = 0;
                image[index++] = 0;
                image[index++] = 255;
                image[index++] = fixedAlpha;

            } else {
                image[index++] = 255;
                image[index++] = 255;
                image[index++] = 255;
                image[index++] = fixedAlpha;
            }
        }
    }
}
void init () {
    glewExperimental = GL_TRUE;
    glewInit (); // init glew, which takes in charges the modern OpenGL calls (v>1.2, shaders, etc)
    glCullFace (GL_BACK);     // Specifies the faces to cull (here the ones pointing away from the camera)
    glEnable (GL_CULL_FACE); // Enables face culling (based on the orientation defined by the CW/CCW enumeration).
    glDepthFunc (GL_LESS); // Specify the depth test for the z-buffer
    glEnable (GL_DEPTH_TEST); // Enable the z-buffer in the rasterization
    glEnableClientState (GL_VERTEX_ARRAY);
    glEnableClientState (GL_NORMAL_ARRAY);
    glEnableClientState (GL_COLOR_ARRAY);
    glEnable (GL_NORMALIZE);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glLineWidth (2.0); // Set the width of edges in GL_LINE polygon mode
    glClearColor (.5f, .5f, .5f, 0.f); // Background color
    colorResponses.resize (tree.positions ().size ());
    camera.resize (DEFAULT_SCREENWIDTH, DEFAULT_SCREENHEIGHT);
    try {
        glProgram = GLProgram::genVFProgram ("Simple GL Program", "src/shaders/shader.vert", "src/shaders/shader.frag"); // Load and compile pair of shaders

    } catch (Exception & e) {
        cerr << e.msg () << endl;
    }

    //Light config
    lightSources.push_back(LightSource(Vec3f(1,0,0),Vec3f(0,1,0),0.6));
    lightSources.push_back(LightSource(Vec3f(-1,0,0),Vec3f(1,0,0),0.6));
    lightSources.push_back(LightSource(Vec3f(0,0,-1),Vec3f(1,1,1),0.5));
    lightSources.push_back(LightSource(Vec3f(0,0,1),Vec3f(1,1,1),0.5));

    //Smoke texture setup

    glGenTextures(1, &smokeTexture);
    glBindTexture(GL_TEXTURE_2D, smokeTexture);
    smokeImage = new unsigned char[4*WIDTH*HEIGHT];
    genCheckerboard(WIDTH, HEIGHT, smokeImage);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, WIDTH, HEIGHT, 0, GL_RGBA, GL_UNSIGNED_BYTE, smokeImage);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);

    camera.getPos(initialCampos);    
    camera.pixelToRay(0,0, frustrumRays[0]);
    camera.pixelToRay(WIDTH-1,0, frustrumRays[1]);
    camera.pixelToRay(WIDTH-1,HEIGHT-1, frustrumRays[2]);
    camera.pixelToRay(0,HEIGHT-1, frustrumRays[3]);
    cout << "Camera: " << initialCampos << endl;
    for(int i = 0; i < 4; i++)
    {
        cout <<"Ray:" << frustrumRays[i] << endl;
    }
    
}

// EXERCISE : the following color response shall be replaced with a proper reflectance evaluation/shadow test/etc.
void updatePerVertexColorResponse () {
  std::vector<Vec3f> positions = tree.positions();
  std::vector<Vec3f> normals = tree.normals();
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
    //colorResponses[i] = Vec3f(128,128,128);
    
        
  }
}

void renderScene () {
    glProgram->use ();

    updatePerVertexColorResponse ();
    glVertexPointer (3, GL_FLOAT, sizeof (Vec3f), (GLvoid*)(&(tree.positions()[0])));
    glNormalPointer (GL_FLOAT, 3*sizeof (float), (GLvoid*)&(tree.normals()[0]));
    glColorPointer (3, GL_FLOAT, sizeof (Vec3f), (GLvoid*)(&(colorResponses[0])));
    glDrawElements (GL_TRIANGLES, 3*tree.triangles().size(), GL_UNSIGNED_INT, (GLvoid*)((&tree.triangles()[0])));
    
    glProgram->stop();
}

void renderSmoke() {
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    gluOrtho2D(0, WIDTH, HEIGHT, 0);
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();  
    glLoadIdentity();
    glDisable(GL_CULL_FACE);
    glClear(GL_DEPTH_BUFFER_BIT);
    glDisable(GL_DEPTH_TEST);

    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, smokeTexture);
    glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
    glBegin(GL_QUADS);
      glTexCoord2f(0, 0); glVertex2f(0., 0.);
      glTexCoord2f(1, 0); glVertex2f(WIDTH, 0);
      glTexCoord2f(1, 1); glVertex2f(WIDTH, HEIGHT);
      glTexCoord2f(0, 1); glVertex2f(0, HEIGHT);
    glEnd();
    glDisable(GL_TEXTURE_2D);
    // Making sure we can render 3d again
    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
    //glPopMatrix();        ----and this?
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_CULL_FACE);
}
void renderInitialCamera(){
    float d = 10.;

    glLineWidth(20.5); 
    glColor3f(1.,1.,1.);
    glBegin(GL_LINES);
    Vec3f rayEnd, ray1, ray2;
    for(int i = 0; i < 4; i++)
    {
        rayEnd = initialCampos + frustrumRays[i] * d;
        glVertex3fv((GLfloat *) &initialCampos);
        glVertex3fv((GLfloat *) &rayEnd);

        ray1 = rayEnd;
        ray2 = initialCampos + frustrumRays[(i+1)%4] * d;
        glVertex3fv((GLfloat *) &ray1);
        glVertex3fv((GLfloat *) &ray2);
    }
    glEnd();
      
}
void reshape(int w, int h) {
    WIDTH = w;
    HEIGHT = h;
    camera.resize (w, h);
}

void display () {
    glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    camera.apply (); 
    renderInitialCamera();
    renderScene ();
    //renderSmoke ();
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
        unsigned int numOfTriangles = tree.triangles ().size ();
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
    init ();
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

