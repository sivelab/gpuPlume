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

  
  //The first three int variables are the number of boxes in
  //the x,y, and z directions respectively.
  //The array of floats is the two coordinate bounds of the cuboid.
  //The first three values is the lower x,y,z coordinate
  //and the last three is the upper x,y,z coordinate
  CollectionBox(int,int,int,float*);
  
  //Checks the location of the particles each pass.
  //If a particle is inside the collection box, it will add
  //one to the particle count for the indexed position.  
  void seeIfInBox(float,float,float);

  //Calculates the moving Average for each position in the 
  //collection box
  void calculateAvg();

  //Writes the moving average values to standard output.
  void outputAvg();

  cell* cBox;

 private:

  int numBox_x; //number of boxes in the x direction
  int numBox_y; //number of boxes in the y direction
  int numBox_z; //number of boxes in the z direction

  //l - lower; u - upper
  float lx;
  float ux;
  float ly;
  float uy;
  float lz;
  float uz;

  //used to keep track of the number of samples taken
  double n;
  //alpha value used in calculating moving average
  double alpha;

};

#endif // __COLLECTIONBOX_H__
