#ifndef __PLUMECONTROL_H__
#define __PLUMECONTROL_H__

//#include <iostream>
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

class PlumeControl{
 public:
  
  PlumeControl(int, int, int);

  void init(); 
  void display();
   
  float time_step; //time step used for the movement of particles
  int twidth,theight;
  int numInRow;

  ParticleControl* pc;
  DisplayControl* dc;
  ParticleEmitter* pe;

  //An array of the collection boxes
  CollectionBox* cBoxes[3];
  //total number of collection boxes being used
  int num_cBoxes;

  std::list<int> indices; 
  std::list<pIndex> indicesInUse;

  GLuint texid[8]; 
  GLenum texType;
  
  FramebufferObject* fbo;
  Renderbuffer* rb;
  GLuint vertex_buffer;

  bool dump_contents;
  bool emit;
  bool show_particle_visuals;
  //bool collectionBoxes;
  bool output_CollectionBox;


 private:
  
  std::list<pIndex>::iterator iter;

  void setupTextures();
  void initFBO();
  void particleReuse();
      
  Timer* display_clock;
  Timer_t display_time[2];
  Timer_t reuse_time[2]; 
  Timer_t cBox_time[2];

  double startCBoxTime;
  double endCBoxTime;
  double totalTime;
  double averagingTime;
  double avgTime;
  bool firstTime;
  bool endCBox;
  char* output_file;

  int frameCount;
  double lifeTime;

  bool odd;
  bool useRealTime;
  bool reuseParticles;

  GLfloat* pos_buffer;

  int testcase;
  GLSLObject emit_shader;
  GLenum int_format; 
  GLenum int_format_init;
  GLint draw_buffer;

  //Total number of particles released in simulation
  double totalNumPar;

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

};

#endif // __PLUMECONTROL_H__
