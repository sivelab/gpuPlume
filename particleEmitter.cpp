#include "particleEmitter.h"
#include <math.h>
#include <iostream>

ParticleEmitter::~ParticleEmitter(){

}
int ParticleEmitter::EmitParticle(FramebufferObject* fbo, bool odd){
  return 0;

}
void ParticleEmitter::Draw(){

  glPointSize(4.0);
  glColor4f(0.0, 0.0, 1.0, 1.0);
  glBegin(GL_POINTS);
  {
    glVertex3f(xpos,ypos,zpos);
  }
  glEnd();

}
void ParticleEmitter::setPosition(float x, float y, float z){
  xpos = x;
  ypos = y;
  zpos = z;

}
void ParticleEmitter::getPosition(float* x, float*y, float*z, float*r){
  *x = xpos;
  *y = ypos;
  *z = zpos;
  *r = 1.0;
}
void ParticleEmitter::getReleasedPosition(float*x,float*y,float*z){
  *x = xpos;
  *y = ypos;
  *z = zpos;
}
void ParticleEmitter::getIndex(int* x, int*y){
  *x = s;
  *y = t;
}

void ParticleEmitter::setParticleReuse(std::list<pIndex>* ind, float time){
  indicesInUse = ind;
 
  lifeTime = time;
  reuse = true;

}
void ParticleEmitter::setNumToEmit(int num){
  //numToEmit = (int)floor(np/tts);
  numToEmit = num;
  
}

//Determines whether or not it is time to emit a particle
bool ParticleEmitter::timeToEmit(float time_step){

  emitTime += time_step*releaseRate;

  if(emitTime >= 1.0){
    //Sets number of particles to emit
    numToEmit = (int)floor(emitTime);
    //Keeps track of the remainder value for emitTime.
    remTime += emitTime - floor(emitTime);

    //When the remainder gets greater than 1, release an extra particle.
    //This is an attempt to be more accurate in releasing particles,
    //instead of just "throwing away" the roundoff values
    if(remTime >= 1.0){
      numToEmit += (int)floor(remTime);
      remTime -= floor(remTime);
    }
    emitTime = 0;
    return true;
  }
  else return false;

}
