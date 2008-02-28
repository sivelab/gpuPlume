#ifndef __VISUALPLANE_H__
#define __VISUALPLANE_H__


#include "particleControl.h"
#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include <list>


typedef struct{
  float x;
  float y;
  float z;
}vec3;


class VisualPlane{

 public:

  VisualPlane(ParticleControl*, float*, float*, float*, float*);
  
  void drawAxisAlignedPlane();
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
  void drawRotationalPlane();
 
  void increaseYaw();
  void increasePitch();
  void increaseRoll();
  void decreaseYaw();
  void decreasePitch();
  void decreaseRoll();

  int plane_layer;
  
  float plane_layer_z,plane_layer_y,plane_layer_x;

  int slider_x;
  bool localValues;
  int visual_field;

  bool move_slider;

  bool rotationPlane;

  //Controls what value the middle color is
  //represented at in the scale used for visualizing
  //the turbulence layers.
  float slider;
  
 private:

  void getLocalTauValues();
  vec3 crossProduct(vec3,vec3);
  float dotProduct(vec3,vec3);
  void Normalize(vec3*);

  void getTextureCoordinates();
  void sortIntersectionPoints();
  void sortTextureCoordinates();
  void calculateNormal();
  bool checkTexCoord();
 
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
  int nxdx,nydy,nzdz;

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
  
  //Point on rotating plane
  vec3 p;

  //normal vector of plane
  vec3 n;

  //offset of plane
  float d;

    
  ///////////////////////
  //Planes of the domain
  ///////////////////////

  //normals
  vec3 n0;
  vec3 n1;
  vec3 n2;
  vec3 n3;
  vec3 n4;
  vec3 n5;
  //offsets
  float d0,d1,d2,d3,d4,d5;

  /////////////////////////

  float eps;

  //Points of rotating plane
  vec3 p0,p1,p2,p3;
  //Texture coordinates for rotating plane
  vec3 t0,t1,t2,t3;

  //List of points found when intersecting planes
  std::list<vec3> pList;
  std::list<vec3> piList;
  std::list<vec3> tList;
  std::list<vec3> tiList;
  std::list<vec3>::iterator listIter;
  std::list<vec3>::iterator texlistIter;
    
  int num_Points;
  int num_Coord;

  float yaw,pitch,roll;


};

#endif //__VISUALPLANE_H__
