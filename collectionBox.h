#ifndef __COLLECTIONBOX_H__
#define __COLLECTIONBOX_H__

#include <math.h>

typedef struct{
  float moving_avg;
  int count;

}cell;

class CollectionBox{
 public:

  CollectionBox(int,int,int,float*);
  
  void seeIfInBox(float,float,float);
  void calculateAvg();

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

  float n;
  float alpha;

};

#endif // __COLLECTIONBOX_H__
