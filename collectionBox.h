#ifndef __COLLECTIONBOX_H__
#define __COLLECTIONBOX_H__
#include <iostream>
#include <math.h>

typedef struct{
  double moving_avg;
  int count;

}cell;

class CollectionBox{
 public:

  CollectionBox(int,int,int,float*);
  
  void seeIfInBox(float,float,float);
  void calculateAvg();
  void outputAvg();

  cell* cBox;

 private:

  int numBox_x;
  int numBox_y;
  int numBox_z;

  float lx;
  float ux;
  float ly;
  float uy;
  float lz;
  float uz;

  double n;
  double alpha;

};

#endif // __COLLECTIONBOX_H__
