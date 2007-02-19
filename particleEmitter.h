#include <GL/glew.h>
#include <GL/glut.h>
#include <iostream>
#include "framebufferObject.h"
#include "GLSL.h"

class ParticleEmitter{

 public:

  ParticleEmitter(float,float,float,int*,int*,std::list<int>*);

  void EmitParticle(FramebufferObject*, GLSLObject, bool);
  void setNumToEmit(int);

 private:

  int numToEmit;

  float xpos,ypos,zpos;

  int twidth,theight;

  std::list<int>* indices;

};
