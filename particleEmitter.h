#ifndef __PARTICLE_EMITTER_H__
#define __PARTICLE_EMITTER_H__
#include <list>
#include <vector>
#include <GL/glew.h>
#include <math.h>
#include "framebufferObject.h"
#include "GLSL.h"

typedef struct{
    int id;
    double time;

}pIndex;


class ParticleEmitter{

 public:

  //The first three values are the x,y,z position
  //The next value is the rate in how many particles to emit per second
  //The two int*'s are the width and height values(which are the two values
  //that determine the number of particles in the system).
  //The list pointer is the list that holds the available indices of particles
  //that can be emitted

  //ParticleEmitter(float,float,float,float,int*,int*,std::list<int>*, GLSLObject*);

  //Emits numToEmit(number of particles to emit) particles by setting values in
  //the particle position textures to xpos,ypos,zpos.
  virtual int EmitParticle(FramebufferObject*, bool);

  //Uses the variables, emitTime and remTime, to emit particles based on the time step
  //and how many particles per second(pps). 
  virtual bool timeToEmit(float);

  //The base function draws the particle emitter as a point
  virtual void Draw();

  virtual void setPosition(float,float,float);

  virtual void getPosition(float*,float*,float*);

  virtual void getReleasedPosition(float*,float*,float*);

  virtual void setParticleReuse(std::list<pIndex>*, float time);

  virtual void setNumToEmit(int);

  virtual void getIndex(int*,int*);

  virtual void setVertices();

  virtual ~ParticleEmitter();

  void setPosTexID(GLuint id0, GLuint id1) { m_posTexID0 = id0; m_posTexID1 = id1; }

  //The position of the particle emitter
  float xpos,ypos,zpos;

  bool releasePerTimeStep;
  bool releaseOne;
  bool releasePerSecond;
  bool instantRelease;
  bool emit;

  bool Punch_Hole;

 protected:

  std::list<float> posCoord;

  GLuint m_posTexID0, m_posTexID1;

  bool reuse;

  //Lifetime of particle
  float lifeTime;

  //Number of particles is twidth*theight
  int twidth,theight;

  //Value used to decide how many particle to emit once the 
  //function, EmitParticle, is called if using particle per second.
  //If particles per time step, numToEmit is set to the number
  //of particles to emit per time step.
  int numToEmit;
  int temp;

  //Index of particle in position texture
  int s;
  int t;
  
  //Release rate of particles:
  //number of particles per second
  float releaseRate;
  
  //These valuse are used in the timeToEmit function,
  //used to determine if it's time to emit a particle or not
  float emitTime, remTime;

  //A list of available particles that have not yet been emitted
  std::list<int>* indices;

  std::list<pIndex>* indicesInUse;

  //A reference to the shader program used to set the initial
  //starting point of the particle
  GLSLObject* shader;

};

#endif // __PARTICLE_EMITTER_H__
