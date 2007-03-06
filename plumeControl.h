#include <GL/glew.h>
#include <GL/glut.h>
#include <iostream>
#include "particleControl.h"
#include "displayControl.h"
#include "particleEmitter.h"
#include "GLSL.h"
#include <list>

class PlumeControl{


 public:
  
  PlumeControl(int, int, int);

  void init();
  void injectParticles(FramebufferObject*, bool);
  void advectParticles(FramebufferObject*, bool);
  void displayVisual(GLuint);
  
 
  float time_step; //time step used for the movement of particles
  int twidth,theight;
  int numInRow;

  ParticleControl* pc;
  DisplayControl* dc;
  ParticleEmitter* pe;
  std::list<int> indices; 

  GLuint texid[8]; 
  GLenum texType;
  

 private:

  void setupTextures();
      
  int testcase;
  GLSLObject emit_shader;
  GLenum int_format; 
  GLenum int_format_init;
  

};
