#ifndef __COLLECTIONBOX_H__
#define __COLLECTIONBOX_H__

class CollectionBox{
 public:

  CollectionBox(int,int,int,float*);
  
  void seeIfInBox(float,float,float);

  int* cBox;

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


};

#endif // __COLLECTIONBOX_H__
