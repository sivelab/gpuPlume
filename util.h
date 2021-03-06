#ifndef __UTIL_H__
#define __UTIL_H__

#include <iostream>
#include <fstream>
#include <sstream>

#include "quicutil/QUBuildings.h"
#include "quicutil/QUMetParams.h"
#include "quicutil/QUSimparams.h"

#include "quicutil/QPBuildout.h"
#include "quicutil/QPParams.h"
#include "quicutil/QPSource.h"

class Util {

 public:

  Util();
  bool readInput(std::string);

  //number of particles
  int twidth,theight;

  //path line dimensions
  int pwidth,pheight;

  float time_step;
  //domain
  int nx,ny,nz;

  //
  float dx,dy,dz;

  //wind field data
  //look in Settings/input.txt file for description of what value to use
  int windFieldData;
  //use real time or fixed time step
  bool useRealTime;
  //output file to store concentration values
  std::string output_file;
  //output id to use when specifying the concentrations in the concentration file
  // will default to "concentration"
  std::string output_id;
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

  // toggle particle reuse
  bool reuse_particles;

  //toggle visualization
  bool show_particle_visuals;
  //toggle collection box visuals
  bool show_collectionBox_visuals;

  //initializes pause mode
  bool pauseMode;

  //Set to do mean velocity calculations
  bool calculateMeanVel;

  //Set to actually perform the colorization step in which the mean velocitys are used to color particles.
  bool updateParticleColors;

  //type of advection ; i.e. determines which shaders and textures to load.
  int advectChoice;

  //total number of particle emitters
  int numOfPE;
  //particle emitter information
  int* petype;
  float* xpos;
  float* ypos;
  float* zpos;
  float* xpos_e;
  float* ypos_e;
  float* zpos_e;
  float* radius;
  float* rate;
  //particle emitter release method
  int releaseType;

  //method for emitting particles
  int emit_method;

  //building paramters
  int* numSides;
  float* xfo;
  float* yfo;
  float* zfo;
  float* ht;
  float* wti;
  float* lti;
  float* gamma;
  int numBuild;
  float calculatedMaxVel;

  bool hasAbsolutePath;
  std::string quicFilesPath;

  //holds the quicPlume Data for the wind field
  //double* u;
  //double* v;
  //double* w;

  //Background color
  float bcolor[3];

  //Number of contour regions
  int num_contour_regions;
  void volumeBox(); //calc Volume of the collecting boxes Balli
  float volume; //volume of the box

  // The mode for the NetworkManager
  int network_mode;

  // The viewing mode for DisplayControl
  int viewing_mode;

  // The view that hsould be used when viewing_mode is set to TREADPORT
  char treadport_view;

  // Value representing wether to use a static or dynamic frustum with the
  // treadport.
  int static_treadport_frustum;

  // fullscreen is an option that specifies if we should run in fullscreen
  // mode or not.
  bool fullscreen;
  
  float sun_azimuth;
  
  float sun_altitude;
  
  bool onlyCalcShadows;

  bool ignoreSignal;
  bool offscreenRender;
  
  // Structures to contain the QU_simparams.inp, QU_buildings.inp, QU_metparams.inp files.
  quSimParams quSimParamData;
  quBuildings quBuildingData;
  quMetParams quMetParamData;

  // A structure to contain the QP_params.inp, QP_source.inp files.
  qpBuildout qpBuildoutData;
  qpParams qpParamData;
  qpSource qpSourceData;

  // problem info
  unsigned int problemID;
  unsigned int problemInstanceID;

  int winWidth, winHeight;

 private:

  void parseLine(char*);
  bool read1Float(char*,std::string,float*);
  bool read6Float(char*,std::string,float*);
  bool read7Float(char*,std::string,float*);
  bool read3Float(char*,std::string,float*);
  bool readSourceInfo(char*,std::string,int&,float*);
  bool readComment(const char*);
  bool read1String(const char*, const char*,std::string*);

  bool isPathAbsolute(const std::string &filename);
  bool isQUICProjFile(std::ifstream& inputStream);
  bool readQUICBaseFiles( std::string& QUICFilesPath );

  //keeps track of number of sources
  int num;
  //keeps track of number of buildings
  int numb;

};

#endif //__UTIL_H__
