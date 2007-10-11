#ifndef __PLUMECONTROL_H__
#define __PLUMECONTROL_H__

#include <iostream>
#include <list>
#include <vector>

#include <GL/glew.h>
#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include "gpuPlume.h"
#include "particleControl.h"
#include "displayControl.h"
#include "particleEmitter.h"
#include "pointEmitter.h"
#include "sphereEmitter.h"
#include "collectionBox.h"
#include "framebufferObject.h"
#include "renderbuffer.h"
#include "GLSL.h"
#include "util.h"
#include "simulation.h"
#include "streamLine.h"
#include "PathLine.h"
#include "Random.h"
#include <string>

class PlumeControl{
 public:

  //PlumeControl(Util*);

  virtual void init(bool); 
  virtual int display();
  virtual void setupTextures();
  virtual void setupEmitters();
  virtual ~PlumeControl();
     
  float time_step; //time step used for the movement of particles
  bool useRealTime; //Set whether to use real time or not

  int twidth,theight; //they are twidth*theight number ofparticles
  int numInRow;
  //Domain 
  int nx;
  int ny;
  int nz;

  StreamLine* stream;
  PathLine* pathLines;

  ParticleControl* pc;
  DisplayControl* dc;
  ParticleEmitter* pe[10];

  //An array of the collection boxes 
  CollectionBox* cBoxes[3];
  //total number of collection boxes being used
  int num_cBoxes;

  std::list<int> indices; 
  std::list<pIndex> indicesInUse;

  GLuint texid[17]; 
  GLenum texType,positions0,positions1,windField,randomValues;
  GLenum prime0, prime1, lambda, tau_dz, duvw_dz, tau;
  GLenum meanVel0,meanVel1,currVel,buildings,cellType;
  GLenum pathTex;

  FramebufferObject* fbo;
  FramebufferObject* fbo2;
  FramebufferObject* pathFbo;
  Renderbuffer* rb;

  GLuint vbo_buffer[2];
  GLenum vertex_buffer;
  GLenum color_buffer;

  bool dump_contents;
  bool print_MeanVel;
  bool createImages;
  bool output_CollectionBox;
  bool osgPlume;
  bool quitSimulation;

  //glut window id
  int winid;

  void getSourceInfo(float*x,float*y,float*z,float*r){
    x = util->xpos;
    y = util->ypos;
    z = util->zpos;
    r = util->radius;   
  }
 
  Simulation* sim;
  float* buildParam;

  bool paused;
  bool inPauseMode;

  //Determines which texture to draw in horizontal layers.
  //Currently the wind field and tau textures are the only options.
  int visualLayer;

  std::string texFileName;

 protected:
  
  GLint currentbuffer,readbuffer;
  std::list<pIndex>::iterator iter;

  virtual void initFBO();
  void particleReuse();
      
  Timer* display_clock;
  Timer_t reuse_time[2]; 
  
  double avgTime;
  bool firstTime;
  bool endCBox;

  int frameCount;
  double lifeTime;

  bool odd;
  bool reuseParticles;

  GLSLObject emit_shader;
  GLenum int_format; 
  GLenum int_format_init;
  GLint draw_buffer;
  GLint read_buffer;

  //Total number of particles released in simulation
  double totalNumPar;

  Util* util;

  GLint* vp;
  GLfloat* mvm;
  GLfloat* pm;

  int maxColorAttachments;

  GLenum currVelBuffer;

};

#endif // __PLUMECONTROL_H__
