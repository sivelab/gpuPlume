#ifndef __IS_H__
#define __IS_H__

#include <GL/glew.h>

#include "util.h"
#include "particleControl.h"
#include "framebufferObject.h"


class IsoSurface{

 public:

  IsoSurface(ParticleControl*);

  void render3DTexture(FramebufferObject*);
  void createIsoSurface();
  void draw();
  void increaseMesh();
  void decreaseMesh();

  float contourValue;
 
  GLuint tex3d[2];

  bool solid;

  bool once;

 private:
  
  void readInTables();

  int nx,ny,nz;
  int nxdx,nydy,nzdz;

  GLenum int_format;

  GLint u_slice, uniform_tau, u_tau3D, u_case, u_edge;
  GLint u_dx,u_dy,u_dz,u_mesh, uniform_cValue;

  GLSLObject geomShader, render3DShader, isoShader;
  
  GLenum texType2;
  
  GLuint case_to_numpoly[1];
  GLuint edge_connect_list[1];
  
  GLuint query;
  int numPrimitives;
  int mesh;

  GLuint iso_buffer[10];
  int buffer_num;

};


#endif
