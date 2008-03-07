#ifndef __CONTOUR_H__
#define __CONTOUR_H__

#include <GL/glew.h>
#include "particleControl.h"
#include <list>

typedef struct{
  float x;
  float y;
}vec2;

class Contour{

 public:
  Contour(ParticleControl*, int);

  void draw();
  void decreaseContourLayer();
  void increaseContourLayer();
  
  void displayContourLayer(ParticleControl*,GLuint,int);

  void switchPlane();

  int tauValue;

 private:
  int nx,ny,nz;
  int nxdx,nydy,nzdz;

  float cell_dx,cell_dy,cell_dz;

  float* cValue;
  int num_cValue;
  int n;

  void setLocalTauValues();
  float* tauLocalMax;
  float* tauLocalMin;
  int max_layer;
  int plane_normal;
  ParticleControl* pc;
  int plane_layer_z,plane_layer_x,plane_layer_y;

  //Stores all the contour values for each height value
  float* contourValues;
  
  GLuint tex_id[3];
  GLuint tex_3D[1];

  float v0,v1,v2,v3;
  vec2 p0,p1,p2,p3;

  float* tau;

  int* numPoints;
  int* numPoints_nx;
  int* numPoints_ny;

  int localPoints;

  std::list<vec2> c1List;
  std::list<vec2>::iterator listIter;

  //VBO's of contours
  GLuint* contourLayer;
  GLuint* contourLayer_nx;
  GLuint* contourLayer_ny;

  //Height layer of Contours to draw
  int layer;

  GLSLObject contour_shader;

  GLint uniform_numContours, uniform_tauTex;
  GLint uniform_tauValue, uniform_contourTex, uniform_height;
  GLint uniform_3Dtau;

  void findContours_Averaging(ParticleControl*);
  void find_Multiple_Contours();
  void find_Multiple_Contours_nx();
  void find_Multiple_Contours_ny();

  void setContourValuesLocally(int);
  void setContourValuesGlobally(ParticleControl*,int);

  void putListinVBO(int);

};

#endif //__CONTOUR_H__
