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

 private:
  int nx,ny,nz;

  float c1;
  float* cValue;
  int num_cValue;

  float v0,v1,v2,v3;
  vec2 p0,p1,p2,p3;

  int numPoints;
  int localPoints;

  std::list<vec2> c1List;
  std::list<vec2>::iterator listIter;

  void findContours_Averaging(ParticleControl*);
  void findContours(ParticleControl*);
  void find_Multiple_Contours(ParticleControl*);

};

#endif //__CONTOUR_H__
