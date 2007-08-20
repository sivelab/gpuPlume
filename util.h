#ifndef __UTIL_H__
#define __UTIL_H__

#include <iostream>
#include <fstream>
#include <sstream>

class Util{

 public:
  
  Util();
  void readInput(std::string);

  //number of particles
  int twidth,theight;

  float time_step;
  //domain
  int nx,ny,nz;
  //wind field data
  //look in Settings/input.txt file for description of what value to use
  int windFieldData;
  //use real time or fixed time step
  bool useRealTime;
  //output file to store concentration values
  std::string output_file;
  //simulation duration in seconds
  double duration;
  //concentration boxes start time
  double startCBoxTime;
  //concentration boxes end time
  double endCBoxTime;
  //concentration averaging time
  double averagingTime;
  //concentration boxes bounds
  float* bounds;
  //number of conc. boxes in x-direction
  int numBox_x;
  //number of conc. boxes in y-direction
  int numBox_y;
  //number of conc. boxes in z-direction
  int numBox_z;
  
  //values for initializing the prime textures
  float ustar,sigU,sigV,sigW;
  
  //toggle visualization
  bool show_particle_visuals;
  //toggle collection box visuals
  bool show_collectionBox_visuals;
  
  //initializes pause mode
  bool pauseMode;

  //Set to do mean velocity calculations
  bool calculateMeanVel;

  //type of advection ; i.e. determines which shaders and textures to load.
  int advectChoice;

  //total number of particle emitters
  int numOfPE;
  //particle emitter information
  float* xpos;
  float* ypos;
  float* zpos;
  float* radius;
  float* rate;
  //particle emitter release method
  int releaseType;
  
  //method for emitting particles
  int emit_method;

  //building paramters
  double* xfo;
  double* yfo;
  double* zfo;
  double* ht;
  double* wti;
  double* lti;
  int numBuild;

  //holds the quicPlume Data for the wind field
  double* u;
  double* v;
  double* w;

  //Background color
  float bcolor[3];

 private:

  void parseLine(char*);
  bool read1Float(char*,std::string,float*);
  bool read6Float(char*,std::string,float*);
  bool read3Float(char*,std::string,float*);
  bool readSourceInfo(char*,std::string,float*);
  bool readComment(const char*);
  bool read1String(const char*,char*,std::string*);

  //keeps track of number of sources
  int num;

};

#endif //__UTIL_H__
