#ifndef __LINE_EMITTER_H__
#define __LINE_EMITTER_H__

#include "particleEmitter.h"

class LineEmitter : public ParticleEmitter{

 public:

  // Arguments to the lineEmitter follow:
  // starting X, starting Y, starting Z, ending X, ending Y, ending Z, 
  // rate, texture width, texture height, indices, emit shader reference, 
  // sig, nxdx, nxdy, nzdz
  LineEmitter(float, float, float, 
	      float, float, float, 
	      float,
	      int,int,std::list<int>*,GLSLObject*,std::vector<float>*,wind*,int,int,int);

  virtual int EmitParticle(bool,GLuint,GLuint,float,GLuint,GLuint);
  virtual void setVertices();
  virtual void Draw();
  virtual void getReleasedPosition(float*,float*,float*);
 
  // we can use the xpos, ypos, and zpos variables from the base class
  // to hold the starting vertex of the line
  float xpos_end, ypos_end, zpos_end;

 protected:
  virtual ~LineEmitter();  


 private:

  float offsetx;
  float offsety;
  float offsetz;

};

#endif //__LINE_EMITTER_H__
