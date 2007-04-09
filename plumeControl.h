#ifndef __PLUMECONTROL_H__
#define __PLUMECONTROL_H__

#include <iostream>
#include <list>

#include <GL/glew.h>
#include <GL/glut.h>

#include "particleControl.h"
#include "displayControl.h"
#include "particleEmitter.h"
#include "pointEmitter.h"
#include "sphereEmitter.h"
#include "collectionBox.h"
#include "gpuPlume.h"
#include "framebufferObject.h"
#include "renderbuffer.h"
#include "GLSL.h"
#include "util.h"
#include "simulation.h"

class PlumeControl{
 public:
  
  PlumeControl(int, int, int);

  void init(bool); 
  void display();
     
  float time_step; //time step used for the movement of particles
  bool useRealTime; //Set whether to use real time or not
  
  int twidth,theight; //they are twidth*theight number ofparticles
  int numInRow;

  ParticleControl* pc;
  DisplayControl* dc;
  ParticleEmitter* pe;

  //An array of the collection boxes
  CollectionBox* cBoxes[3];
  //total number of collection boxes being used
  int num_cBoxes;

  //testcase determines which data set to use for the windfield.
  //The value t is currently passed in from gpuPlume.  When it
  //equals 3, it runs the quicplume data set.  When it equals
  //4 it runs the uniform u-direction windfield.  
  int testcase;

  std::list<int> indices; 
  std::list<pIndex> indicesInUse;

  GLuint texid[8]; 
  GLenum texType,positions0,positions1,windField,randomValues;
  
  FramebufferObject* fbo;
  Renderbuffer* rb;
  GLuint vertex_buffer;

  bool dump_contents;
  bool emit;
  bool show_particle_visuals;
  bool output_CollectionBox;
  bool osgPlume;
  bool quitSimulation;
 
  //output file of collection box
  std::string output_file;
  
  //simulation duration in seconds
  double duration;

  //glut window id
  int winid;

 private:
  
  Simulation* sim;
  
  std::list<pIndex>::iterator iter;

  void setupTextures();
  void initFBO();
  void particleReuse();
  void readInputFile();
      
  Timer* display_clock;
  Timer_t reuse_time[2]; 
  
  double startCBoxTime;
  double endCBoxTime;
  double averagingTime;
  double avgTime;
  bool firstTime;
  bool endCBox;

  int frameCount;
  double lifeTime;

  bool odd;
  bool reuseParticles;
  bool releasePerTimeStep;

  GLfloat* pos_buffer;
  GLSLObject emit_shader;
  GLenum int_format; 
  GLenum int_format_init;
  GLint draw_buffer;

  //Total number of particles released in simulation
  double totalNumPar;

  Util* utility;

  //QUIC-PLUME References
  int nx;
  int ny;
  int nz;
  double* u;
  double* v;
  double* w;

  int numBuild;
  double* xfo;
  double* yfo;
  double* zfo;
  double* ht;
  double* wti;
  double* lti;


  GLint* vp;
  GLfloat* mvm;
  GLfloat* pm;

};

#endif // __PLUMECONTROL_H__
