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
  Contour(ParticleControl*);

  void draw();
  void decreaseContourLayer();
  void increaseContourLayer();
  
  void displayContourLayer(ParticleControl*,GLuint,int);

  int tauValue;

 private:
  int nx,ny,nz;

  float* cValue;
  int num_cValue;
  int n;

  //Stores all the contour values for each height value
  float* contourValues;
  
  GLuint tex_id[1];

  float v0,v1,v2,v3;
  vec2 p0,p1,p2,p3;

  float* tau;

  int* numPoints;
  int localPoints;

  std::list<vec2> c1List;
  std::list<vec2>::iterator listIter;

  //VBO's of contours
  GLuint* contourLayer;

  //Height layer of Contours to draw
  int layer;

  GLSLObject contour_shader;

  GLint uniform_numContours, uniform_tauTex;
  GLint uniform_tauValue, uniform_contourTex, uniform_height;
  
  void findContours_Averaging(ParticleControl*);
  void find_Multiple_Contours(ParticleControl*);

  void setContourValuesLocally(ParticleControl*,int);
  void setContourValuesGlobally(ParticleControl*,int);

  void putListinVBO(int);

};

#endif //__CONTOUR_H__
