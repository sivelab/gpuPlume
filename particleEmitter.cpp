#include "particleEmitter.h"

ParticleEmitter::~ParticleEmitter(){

}
void ParticleEmitter::EmitParticle(FramebufferObject* fbo, bool odd){

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
void ParticleEmitter::setParticleReuse(std::list<pIndex>* ind, float time){
  indicesInUse = ind;
 
  lifeTime = time;
  reuse = true;

}

//Determines whether or not it is time to emit a particle
bool ParticleEmitter::timeToEmit(float time_step){
  emitTime += time_step*pps;

  if(emitTime >= 1.0){
    remTime += emitTime - floor(emitTime);

    if(remTime >= 1.0){
      numToEmit += (int)floor(remTime);
      remTime -= floor(remTime);
    }
    emitTime = 0;
    return true;
  }
  else return false;

}
