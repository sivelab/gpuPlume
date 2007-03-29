#ifndef __POINT_EMITTER_H__
#define __POINT_EMITTER_H__

#include "particleEmitter.h"

class PointEmitter : public ParticleEmitter{

 public:

  PointEmitter(float,float,float,float,int*,int*,std::list<int>*, GLSLObject*);
  virtual int EmitParticle(FramebufferObject*, bool);
  //virtual void Draw();

 protected:
  virtual ~PointEmitter();

};

#endif //__POINT_EMITTER_H__
