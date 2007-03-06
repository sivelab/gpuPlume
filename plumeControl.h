#ifndef __PLUMECONTROL_H__
#define __PLUMECONTROL_H__

#include <iostream>
#include <list>

#include <GL/glew.h>
#include <GL/glut.h>
#include "particleControl.h"
#include "displayControl.h"
#include "particleEmitter.h"
#include "framebufferObject.h"
#include "renderbuffer.h"
#include "GLSL.h"

class PlumeControl{


 public:
  
  PlumeControl(int, int, int);

  void init();
 
  void display();
  //void injectParticles(FramebufferObject*, bool);
  //void advectParticles(FramebufferObject*, bool);
  //void displayVisual(GLuint);
  
 
  float time_step; //time step used for the movement of particles
  int twidth,theight;
  int numInRow;

  ParticleControl* pc;
  DisplayControl* dc;
  ParticleEmitter* pe;
  std::list<int> indices; 

  GLuint texid[8]; 
  GLenum texType;
  
  FramebufferObject* fbo;
  Renderbuffer* rb;
  GLuint vertex_buffer;

  bool dump_contents;
  bool emit;
  bool show_particle_visuals;

 
 

 private:

  void setupTextures();
  void initFBO();
      
  Timer* display_clock;
  Timer_t display_time[2];

  bool odd;

  int testcase;
  GLSLObject emit_shader;
  GLenum int_format; 
  GLenum int_format_init;

  GLint draw_buffer;
};

#endif // __PLUMECONTROL_H__
