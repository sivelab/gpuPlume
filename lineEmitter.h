#ifndef __LINE_EMITTER_H__
#define __LINE_EMITTER_H__

#include "particleEmitter.h"

class LineEmitter : public ParticleEmitter{

 public:

  LineEmitter(float,float,float,float,float,int,int,std::list<int>*,GLSLObject*);

  virtual int EmitParticle(bool,GLuint,GLuint,float);
  virtual void setVertices();
  virtual void Draw();
  virtual void getReleasedPosition(float*,float*,float*);
 
 protected:
  virtual ~LineEmitter();  

 private:
  float radius;

  float offsetx;
  float offsety;
  float offsetz;

};

#endif //__LINE_EMITTER_H__
