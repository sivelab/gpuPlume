#ifndef __COLLECTIONBOX_H__
#define __COLLECTIONBOX_H__
#include <iostream>
#include <math.h>
#include <fstream>

/*typedef struct{
  double concentration;
  int count;

  }cell;*/

class CollectionBox{
 public:

  
  //The first three int variables are the number of boxes in
  //the x,y, and z directions respectively.
  //The array of floats is the two coordinate bounds of the cuboid.
  //The first three values is the lower x,y,z coordinate
  //and the last three is the upper x,y,z coordinate
  //The last float is the concentration averaging time
  CollectionBox(int,int,int,float*,double);
  
  //Checks the location of the particles each pass. 
  //void seeIfInBox(float,float,float);

  //Calculates the Concentration for each position in the 
  //collection box
  void calculateConc(float,float,float,float,double);
  
  //Set concentration values to zero
  void clear();

  //Writes the concentration values to specified file
  void outputConc(std::string,double);

  //cell* cBox;
  double* cBox;

 private:

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


};

#endif // __COLLECTIONBOX_H__
