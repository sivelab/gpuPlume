#ifndef __COLLECTIONBOX_H__
#define __COLLECTIONBOX_H__
#include <iostream>
#include <math.h>
#include <fstream>
#include <GL/glew.h>
#include "simulation.h"
#include "util.h"

typedef struct{
  float x;
  float y;
  float z;
}BoxPos;

typedef struct{
  float d;
  int idx;
}BoxDis;

class CollectionBox{
 public:
  
  CollectionBox(Util*);

  bool findConc(Simulation*,bool*,bool);

  void calcSimpleConc(float,float,float);

  //Calculates the Concentration for each position in the 
  //collection box
  void calculateConc(float,float,float,float,double);
  
  //Set concentration values to zero
  void clear();
  void draw(double);
  void sort(float,float,float);

  //Writes the concentration values to specified file
  void outputConc(std::string,double,double);

  double* cBox; 
  BoxPos* cPos;
  BoxDis* cDis;

  float dx;
  float dy;
  float dz;
 private:

  void calcBoxPositions();

  //std::ofstream output;
  bool alreadyOpen;

  int numBox_x; //number of boxes in the x direction
  int numBox_y; //number of boxes in the y direction
  int numBox_z; //number of boxes in the z direction

  //l = lower; u = upper
  float lx;
  float ux;
  float ly;
  float uy;
  float lz;
  float uz;

  double volume;
  double constant;
  double TotRel;
  double concAvgTime;

  double endCBoxTime,startCBoxTime,avgTime;

  int twidth,theight;
  GLfloat* pos_buffer;

};

#endif // __COLLECTIONBOX_H__
