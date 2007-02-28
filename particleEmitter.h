#include <GL/glew.h>
#include <GL/glut.h>
#include <iostream>
#include "framebufferObject.h"
#include "GLSL.h"

class ParticleEmitter{

 public:

  //The first three values are the x,y,z position
  //The next value is the rate in how many particles to emit per second
  //The two int*'s are the width and height values(which are the two values
  //that determine the number of particles in the system).
  //The list pointer is the list that holds the available indices of particles
  //that can be emitted.
  ParticleEmitter(float,float,float,float,int*,int*,std::list<int>*, GLSLObject*);

  //Emits numToEmit(number of particles to emit) particles by setting values in
  //the particle position textures to xpos,ypos,zpos.
  void EmitParticle(FramebufferObject*, bool);

  //Uses the variables, emitTime and remTime, to emit particles based on the time step
  //and how many particles per second(pps). 
  bool timeToEmit(float);

 private:

  //The position of the particle emitter
  float xpos,ypos,zpos;

  //Number of particles is twidth*theight
  int twidth,theight;

  //Value used to decide how many particle to emit once the function, EmitParticle,
  //is called.  
  int numToEmit;
  
  //Number of particles to emit per second
  float pps;

  float emitTime, remTime;

  std::list<int>* indices;

  GLSLObject* shader;

};
