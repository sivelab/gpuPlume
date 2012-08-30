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
#include "lineEmitter.h"
#include "planeEmitter.h"
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
#include "Contour.h"
#include "VisualPlane.h"
#include "IsoSurface.h"
#include <string>


class PlumeControl{
 public:
 
  //PlumeControl(Util*);

  virtual void init(bool); 
  virtual int display(long int );  ///value being passed over to print within the concentration file 
  virtual void setupTextures();
  virtual void setupEmitters();
  virtual void swapPauseMode();
  virtual ~PlumeControl();
  
  // writeShadowMapToFile will write the generated
  // shadow map (calculated at the beginning of the
  // run) to a file. Currently, the file is a text
  // file and the file name is static.

  // Note, this method should really be in MulitpleBuildings
  // but since we have a virtual plume controle it is here.
  // This means that the other models should implement this
  // method at some point.
  virtual void writeShadowMapToFile();
  
  float time_step; //time step used for the movement of particles
  bool useRealTime; //Set whether to use real time or not

  int twidth,theight; //they are twidth*theight number ofparticles
  int numInRow;
  //Domain 
  int nx;
  int ny;
  int nz;

  int nxdx,nydy,nzdz;

  StreamLine* stream;
  PathLine* pathLines;
  Contour* contours;

  VisualPlane* planeVisual;
  IsoSurface* isoSurface;

  ParticleControl* pc;
  DisplayControl* dc;
  ParticleEmitter* pe[10];

  //An array of the collection boxes 
  CollectionBox* cBoxes[3];
  //total number of collection boxes being used
  int num_cBoxes;

  std::list<int> indices; 
  std::list<pIndex> indicesInUse;

  GLuint texid[19]; 
  GLenum texType,positions0,positions1,windField,randomValues;
    GLenum prime0, prime1, lambda, tau_dz, duvw_dz, tau,dxyz_wall;
  GLenum meanVel0,meanVel1,currDirection,buildings,cellType;
  GLenum pathTex, advect_terms;

  FramebufferObject* fbo;
  FramebufferObject* fbo2;
  FramebufferObject* pathFbo;
  FramebufferObject* isoFbo;
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
  bool drawIsoSurface;
  bool color_by_advect_terms;

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

  std::string texFileName;
  int oneTime;

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
  bool continuousParticleFlow;

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

  GLenum particleColorBuffer;

  std::vector<float> random_values;

};

#endif // __PLUMECONTROL_H__
