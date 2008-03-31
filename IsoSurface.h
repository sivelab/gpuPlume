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
  void createIsoSurface(GLuint);
  void draw(GLuint);
  
  GLuint tex3d[2];
 private:
  
  void readInTables();

  int nx,ny,nz;
  int nxdx,nydy,nzdz;

  GLenum int_format;

  GLint u_slice, uniform_tau, u_tau3D, u_case, u_edge;
  GLint u_dx,u_dy,u_dz,u_mesh;

  GLSLObject geomShader, render3DShader, isoShader;
  
  GLenum texType2;
  
  GLuint case_to_numpoly[1];
  GLuint edge_connect_list[1];
  
  GLuint query;
  int numPrimitives;
  int mesh;

};


#endif
