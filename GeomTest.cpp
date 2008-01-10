#include <math.h>
#include "GeomTest.h"
#include "glErrorUtil.h"


#ifdef WIN32
#include <windows.h>
#include <stdio.h>
#include <conio.h>
#include <tchar.h>

#endif

GeomTest::GeomTest(Util* u){
  
  util = u;

  //from util
  twidth = util->twidth;
  theight = util->theight;
  nx = util->nx;
  ny = util->ny;
  nz = util->nz;
  time_step = util->time_step;

  odd = true;
  dump_contents = true;
  //paused = false;

}
GeomTest::~GeomTest(){}

void GeomTest::init(bool OSG){

  //pc = new ParticleControl(texType, twidth,theight,nx,ny,nz);

  
  int num_vertices;
  if(GL_EXT_geometry_shader4)
    {
      std::cout << "Ready for geom shader!" << std::endl;
      //int numV = 0;
      //glGetIntegerv(GL_MAX_GEOMETRY_OUTPUT_VERTICES_EXT, &numV);
      
      glGetIntegerv(GL_MAX_GEOMETRY_OUTPUT_VERTICES_EXT, &num_vertices);
      std::cout << num_vertices << " number of output vertices." << std::endl;
    }
  geomShader.addShader("Shaders/streamLines_vp.glsl", GLSLObject::VERTEX_SHADER);
  geomShader.addShader("Shaders/streamLines_gp.glsl", GLSLObject::GEOMETRY_SHADER);
  geomShader.addShader("Shaders/streamLines_fp.glsl", GLSLObject::FRAGMENT_SHADER);
  geomShader.setInputandOutput(GL_POINTS,GL_TRIANGLE_STRIP,num_vertices);
  geomShader.createProgram();
  


}

int GeomTest::display(){
    
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  glViewport(0, 0, glutGet(GLUT_WINDOW_WIDTH), glutGet(GLUT_WINDOW_HEIGHT));
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(60.0, glutGet(GLUT_WINDOW_WIDTH)/float(glutGet(GLUT_WINDOW_HEIGHT)), 1.0, 250.0);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  glClearColor(0.0,0.0,0.0,1.0);
	
  gluLookAt(10,0,0,0,0,0, 0, 0, 1 );
  
  geomShader.activate();
  glPointSize(4.0);
  glBegin(GL_POINTS);

  glColor4f(1.0,0.0,0.0,1.0);
  glVertex3f(0.0,0.0,0.0);
  glVertex3f(0.0,1.0,0.0);

  glEnd();

  geomShader.deactivate();

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluOrtho2D(-1, 1, -1, 1);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
     
  //glDisable(texType);
  CheckErrorsGL("END : visualization");


  glutSwapBuffers();
  return 1;

}
