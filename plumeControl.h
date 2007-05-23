#ifndef __PLUMECONTROL_H__
#define __PLUMECONTROL_H__

#include <iostream>
#include <list>
#include <vector>

#include <GL/glew.h>
#include <GL/glut.h>

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
#include "Random.h"

class PlumeControl{
 public:
  
  PlumeControl();

  void init(bool); 
  int display();
     
  float time_step; //time step used for the movement of particles
  bool useRealTime; //Set whether to use real time or not
  
  int twidth,theight; //they are twidth*theight number ofparticles
  int numInRow;
  //Domain 
  int nx;
  int ny;
  int nz;

  StreamLine* stream;
  ParticleControl* pc;
  DisplayControl* dc;
  ParticleEmitter* pe[10];
  int numOfPE;

  //An array of the collection boxes 
  CollectionBox* cBoxes[3];
  //total number of collection boxes being used
  int num_cBoxes;
  double startCBoxTime;
  double endCBoxTime;
  double averagingTime;
  float* bounds;
  int numBox_x, numBox_y, numBox_z;

  //testcase determines which data set to use for the windfield.
  //The value t is currently passed in from gpuPlume.  When it
  //equals 3, it runs the quicplume data set.  When it equals
  //4 it runs the uniform u-direction windfield.  
  int testcase;

  std::list<int> indices; 
  std::list<pIndex> indicesInUse;

  GLuint texid[8]; 
  GLenum texType,positions0,positions1,windField,randomValues;
  GLenum prime0, prime1, lambda;
  
  FramebufferObject* fbo;
  Renderbuffer* rb;
  GLuint vertex_buffer;

  bool dump_contents;
  bool show_particle_visuals;
  bool output_CollectionBox;
  bool osgPlume;
  bool quitSimulation;
  bool show_collectionBox_visuals;
 
  //output file of collection box
  std::string output_file;
  
  //simulation duration in seconds
  double duration;

  //glut window id
  int winid;

  float ustar,sigU,sigV,sigW;
  
  //ParicleEmitter information (source)
  //float xpos,ypos,zpos,radius;
  float* xpos;
  float* ypos;
  float* zpos;
  float* radius;

  int releaseType;

  void getSourceInfo(float*x,float*y,float*z,float*r){
    x = xpos;
    y = ypos;
    z = zpos;
    r = radius;   
  }
 
  Simulation* sim;

 private:
  
  
  std::list<pIndex>::iterator iter;

  void setupTextures();
  void setupEmitters();
  void initFBO();
  void particleReuse();
  //void readInputFile();
      
  Timer* display_clock;
  Timer_t reuse_time[2]; 
  
  double avgTime;
  bool firstTime;
  bool endCBox;

  int frameCount;
  double lifeTime;

  bool odd;
  bool reuseParticles;

  GLfloat* pos_buffer;
  GLSLObject emit_shader;
  GLenum int_format; 
  GLenum int_format_init;
  GLint draw_buffer;

  //Total number of particles released in simulation
  double totalNumPar;

  Util* utility;

  //QUIC-PLUME References
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
