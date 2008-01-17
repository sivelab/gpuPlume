#ifndef __VISUALPLANE_H__
#define __VISUALPLANE_H__


#include "particleControl.h"
#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif


typedef struct{
  float x;
  float y;
  float z;
}vec3;


class VisualPlane{

 public:

  VisualPlane(Matrix*,int,int,int, float*,float*,float*,float*);
  
  void drawPlane();
  void increasePlaneLayer();
  void decreasePlaneLayer();
  bool clickedSlider(int,int);
  void moveSlider(int);
  void moveSliderDown();
  void moveSliderUp();
  void clickedRangeButton(int,int);
  void drawScale();
  void switchPlane();
  void getIntersectionPoints();

  int plane_layer;
  int slider_x;
  bool localValues;
  int visual_field;

  bool move_slider;

  //Controls what value the middle color is
  //represented at in the scale used for visualizing
  //the turbulence layers.
  float slider;

 private:

  void getLocalTauValues();

  float* TauMax;
  float* TauMin;
  float* tauLocalMax;
  float* tauLocalMin;
  float tauMin, tauMax;
  char* Taus[4];

  float scale_xstart,scale_xend;
  float scale_ystart,scale_yend;

  int* tauPos_x;

  int rangeButton_x,rangeButton_y; 

  GLuint tex_id[1];
  GLenum texType, format;

  int nx,ny,nz;

  GLint u_tauTex;
  GLint u_max11,u_max22,u_max33,u_max13;
  GLint u_min11,u_min22,u_min33,u_min13;
  GLint u_controlTau, u_sliderTurb;
  GLint uniform_xmax, uniform_xmin, uniform_tauMin, uniform_tauMax;
  GLint uniform_sliderScale;

  GLSLObject plane_shader, scale_shader;

  Matrix* tau;

  int plane_normal;
  int max_layer;
  
  //Angle for rotating plane
  float thetaX;
  
  //normal vector of plane
  vec3 n;
  //Point on edge of rotating plane
  vec3 r1;
  //Second point used to find intersection point on domain
  vec3 r2;
  //Point on domain plane used to find intersection point
  vec3 V;
  //Normal of the domain plane
  vec3 N;
   
  //Points of rotating plane
  vec3 p0,p1,p2,p3;

};

#endif //__VISUALPLANE_H__
