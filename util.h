#ifndef __UTIL_H__
#define __UTIL_H__

#include <iostream>
#include <fstream>
#include <sstream>

class qpParams
{
public:
  enum SourceType {
    BASIC = 1,
    DENSEGAS = 2,
    DISTPARTSIZE = 3,
    EXPLOSIVE = 4,
    ERADSOURCE = 5,
    BIOSLURRY = 6,
    TWOPHASE = 7,
    EXFILTRATION = 8    
  };
  
  SourceType sourceFlag;  
  short isiteflag;   // !normal QUIC (isitefl=0) or sensor siting (=1) mode
  bool iindoorflag; // !indoor calculations turned off (=0) or turned on (=1)
  short inextgridflag;      // !1 - inner grid, 2 - outer grid
  float westernEdge;  // !Location of western edge of inner grid relative to outer grid (m)
  float southernEdge; // !Location of southern edge of inner relative to outer grid (m)
  float z0;  // wallSurfRoughness;    // !Wall Surface Roughness Length (m)
  float rcl;  // !Reciprocal Monin-Obukhov length(1/m)
  float boundaryLayerHeight;  // !Boundary Layer height (m)
  bool nonLocalMixing;        // !use 1 to enable non-local mixing
  int numParticles;          // !number of particles released over entire simulation
  short particleDistFlag;     // !Number of particle distribution flag (1 = by mass, 2 = by source)
  bool particleSplitFlag;     // !Particle splitting flag
  bool particleRecyclingFlag; // !Particle recycling flag
  int partNumIncreaseFactor;  // !Total particle number increase factor
  short numParticleSplit;     // !Number of particles a particle is split into
  double partSplittingDosage; // !Particle splitting target dose (gs/m^3)
  float taylorMicroscaleMin;  // !Enable Taylor microscale lower limit to sub-time steps
  int randomNumberSeed;  // !Random number seed
  double timeStep;       // !time step (s)
  double duration;       // !duration (s)
  double concAvgTime;   // !concentration averaging time (s)
  double concStartTime; // !starting time for concentration averaging (s)
  double partOutputPeriod; // !particle output period (s)
  float nbx;  // !in x direction, # of collecting boxes (concentration grid cells) 
  float nby;  // !in y direction, # of collecting boxes (concentration grid cells) 
  float nbz;  // !in z direction, # of collecting boxes (concentration grid cells) 
  float xbl;  // !lower limits for collecting boxes in x in meters
  float xbu;  // !upper limits for collecting boxes in x direction in meters
  float ybl;  // !lower limits for collecting boxes in y in meters
  float ybu;  // !upper limits for collecting boxes in y direction in meters
  float zbl;  // !lower limits for collecting boxes in z in meters
  float zbu;  // !upper limits for collecting boxes in z direction in meters
};

class Util{

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
  
  // A structure of parameters from the QP_params.inp file.
  qpParams qpParamData;

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
