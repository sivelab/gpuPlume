#ifndef __POINT_EMITTER_H__
#define __POINT_EMITTER_H__

#include "particleEmitter.h"

class PointEmitter : public ParticleEmitter{

 public:

  PointEmitter(float,float,float,float,int,int,std::list<int>*,GLSLObject*,std::vector<float>*,wind*,
	       int,int,int);
  virtual int EmitParticle(bool,GLuint,GLuint, float, GLuint,GLuint);
  virtual void setVertices();
  //virtual void Draw();

 protected:
  virtual ~PointEmitter();

};

#endif //__POINT_EMITTER_H__
