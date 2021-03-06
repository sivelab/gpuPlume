#ifndef __PARTICLE_EMITTER_H__
#define __PARTICLE_EMITTER_H__
#include <list>
#include <vector>
#include <GL/glew.h>
#include <math.h>
#include "framebufferObject.h"
#include "GLSL.h"
#include "particleControl.h"

struct pIndex { 
    int id;
    double time;
    
    pIndex(const pIndex &cpy) {
      id = cpy.id;
      time = cpy.time;
    }
    
    pIndex() {
      id = 0;
      time = 0.0;
    }
};

enum particleReleaseType{perTimeStep,perSecond,onePerKeyPress,instantaneous};

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
  virtual int EmitParticle(bool, GLuint, GLuint, float,GLuint,GLuint);

  //Uses the variables, emitTime and remTime, to emit particles based on the time step
  //and how many particles per second(pps). 
  virtual bool timeToEmit(float);

  //The base function draws the particle emitter as a point
  virtual void Draw();

  virtual void setPosition(float,float,float);

  virtual void getPosition(float*,float*,float*);

  virtual void getReleasedPosition(float*,float*,float*);

  virtual void setParticleReuse(std::list<pIndex>*, float time);

  virtual void setContinuousParticleFlow(std::list<pIndex>*, float time);

  virtual void setNumToEmit(float);

  virtual void setEmitAmount(float);
  
  virtual int addRemainder();

  virtual void getIndex(int*,int*);

  virtual void setVertices();

  virtual ~ParticleEmitter();

  //void setPosTexID(GLuint id0, GLuint id1) { m_posTexID0 = id0; m_posTexID1 = id1; }

  //The position of the particle emitter
  float xpos,ypos,zpos;

  bool emit;
  bool Punch_Hole;
  particleReleaseType releaseType;
  
 protected:

  std::list<float> posCoord;

  //random values 
  std::vector<float>* random_values; 

  //Balli: included rotation parameter
  //[didn't know how to achieve this using "references", may be faster if we pass a "const" pointer instead of passing whole vector]
  std::vector<float> alph1ij,alph2ij,alph3ij,bet1ij,bet2ij,bet3ij,gam1ij,gam2ij,gam3ij;

  //Sigma values for wind field from particleControl
  wind* sigma;

  int nxdx,nydy,nzdz;
  unsigned int curr;

  GLuint m_posTexID0, m_posTexID1;

  bool reuse;
  bool continuosParticles;

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

  float emitAmount;
  float remAmount;

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
