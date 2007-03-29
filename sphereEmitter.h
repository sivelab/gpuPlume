#ifndef __SHPERE_EMITTER_H__
#define __SPHERE_EMITTER_H__

#include "particleEmitter.h"
#include <GL/glut.h>

class SphereEmitter : public ParticleEmitter{

 public:

  SphereEmitter(float,float,float,float,float,int*,int*,std::list<int>*, GLSLObject*);

  virtual void EmitParticle(FramebufferObject*, bool);
  virtual void Draw();
 
 protected:
  virtual ~SphereEmitter();  

 private:
  float radius;

  float offsetx;
  float offsety;
  float offsetz;

};

#endif //__SPHERE_EMITTER_H__