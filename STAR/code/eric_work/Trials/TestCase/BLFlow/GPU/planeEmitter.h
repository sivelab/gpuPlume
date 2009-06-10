#ifndef __PLANE_EMITTER_H__
#define __PLANE_EMITTER_H__

#include "particleEmitter.h"

class PlaneEmitter : public ParticleEmitter{

 public:

  // Arguments to the planeEmitter follow:
  // origin X, origin Y, origin Z, width, height, 
  // rate, texture width, texture height, indices, emit shader reference, 
  // sig, nxdx, nxdy, nzdz
  PlaneEmitter(float, float, float, 
	       float, float,  
	       float,
	       int,int,std::list<int>*,GLSLObject*,std::vector<float>*,wind*,int,int,int);

  virtual int EmitParticle(bool,GLuint,GLuint,float,GLuint,GLuint);
  virtual void setVertices();
  virtual void Draw();
  virtual void getReleasedPosition(float*,float*,float*);
 
 protected:
  virtual ~PlaneEmitter();  

 private:

  // The origin of the plane is held in the xpos, ypos, and zpos
  // variables from the base class.
  float plane_width, plane_height;
  float offsetx, offsety, offsetz;
};

#endif //__PLANE_EMITTER_H__
