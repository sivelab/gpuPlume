#include "util.h"
#include <stdio.h>
#include "gpuPlume.h"

// //////////////////////////////////////
// BEGIN -----> QUIC PLUME FORTRAN REFERENCES
// //////////////////////////////////////

#if 0

extern "C"
{
  void readfiles_();
}

// Domain size stored in nx, ny, and nz
extern "C" int __datamodule__nx;
extern "C" int __datamodule__ny;
extern "C" int __datamodule__nz;

extern "C" double __datamodule__dx;
extern "C" double __datamodule__dy;
extern "C" double __datamodule__dz;

// UVW contains the wind field
extern "C" double* __datamodule__u;
extern "C" double* __datamodule__v;
extern "C" double* __datamodule__w;

extern "C" int __datamodule__inumbuild;   // integer number of buildings
extern "C" double* __datamodule__xfo;
extern "C" double* __datamodule__yfo;
extern "C" double* __datamodule__zfo; 
extern "C" double* __datamodule__ht;
extern "C" double* __datamodule__wti;
extern "C" double* __datamodule__lti; 

#endif
// //////////////////////////////////////
// END ----> QUIC PLUME FORTRAN REFERENCES
// //////////////////////////////////////

Util::Util(){
  num = 0;
  numb = 0;
  bounds = new float[6];
}

void Util::readInput(std::string file){

  std::ifstream in;
  in.open(file.c_str(),std::ios::in);
    	
  if(in == NULL) std::cout << "input file didn't open" << std::endl;
    
  char line[1024];

  while(  !in.eof() )
  {  
	 in.getline(line, 1024);
	 if( line[ strlen(line)] == '\n' ){
		   line[ strlen(line)] = '\0';
	  }

       parseLine(line);
  }
    
  in.close();

#ifdef USE_PLUME_DATA
 // Call the PLUME code to read in the data files.
  std::cout << "Reading data using PLUME code..." << std::endl;
  readfiles_();

  //QuicPlume data for the domain
  nx = __datamodule__nx; //domain in the x direction
  ny = __datamodule__ny; //domain in the y direction
  nz = __datamodule__nz; //domain in the z direction

  //nx = (__datamodule__nx - 1) * __datamodule__dx; //domain in the x direction
  //ny = (__datamodule__nz - 1) * __datamodule__dz; //domain in the y direction
  //nz = (__datamodule__ny - 1) * __datamodule__dy; //domain in the z direction

  std::cout << "QUIC PLUME domain size: " << nx << " (in X) by " 
	    << ny << " (in Y) by " << nz << " (in Z)" << std::endl;

  //QuicPlume data for the windfield
  u = __datamodule__u;
  v = __datamodule__v;
  w = __datamodule__w;

  //QuicPlume data for the buildings
  numBuild = __datamodule__inumbuild;
  xfo = __datamodule__xfo;
  yfo = __datamodule__yfo;
  zfo = __datamodule__zfo;
  ht = __datamodule__ht;
  wti = __datamodule__wti;
  lti = __datamodule__lti;
  
#else
  //nx = 60;
  //ny = 60;//140;
  //nz = 20;
  u=0;
  v=0;
  w=0;
#endif

}
void Util::parseLine(char* line){
  
  float f1;
  float* b = new float[6];
  float* s = new float[5];
  float* c = new float[3];
  std::string s1;

  if(readComment(line))
    return;

  if(read1Float(line, "twidth", &f1)){
	 twidth = (int)f1;
  }
  if(read1Float(line, "theight", &f1)){
	 theight = (int)f1;
  }
  if(read1Float(line, "time_step", &f1)){
    time_step = f1;
  }
  if(read1Float(line, "nx", &f1)){
    nx = (int)f1;
  }
  if(read1Float(line, "ny", &f1)){
    ny = (int)f1;
  }
  if(read1Float(line, "nz", &f1)){
    nz = (int)f1;
  }
  if(read1Float(line, "windFieldData", &f1)){
    windFieldData = (int)f1;
  }
  if(read1Float(line, "useRealTime", &f1)){
    if(f1 == 0)
      useRealTime = false;
    else
      useRealTime = true;
  }
  if(read1String(line, "output_file", &s1)){
    output_file = s1;
  }
  if(read1Float(line, "duration", &f1)){
    duration = (double)f1;
  }
  if(read1Float(line, "startCBoxTime", &f1)){
    startCBoxTime = (double)f1;
  }
  if(read1Float(line, "endCBoxTime", &f1)){
    endCBoxTime = (double)f1;
  }
  if(read1Float(line, "averagingTime", &f1)){
    averagingTime = (double)f1;
  }
  if(read6Float(line, "bounds", b)){
    bounds = b;
  }
  if(read1Float(line, "numBox_x", &f1)){
    numBox_x = (int)f1;
  }
  if(read1Float(line, "numBox_y", &f1)){
    numBox_y = (int)f1;
  }
  if(read1Float(line, "numBox_z", &f1)){
    numBox_z = (int)f1;
  }
  if(read1Float(line, "ustar", &f1)){
    ustar = f1;
    sigU = 2.0*f1;
    sigV = 2.0*f1;
    sigW = 1.3*f1;
  }
  if(read1Float(line, "show_particle_visuals", &f1)){
    if(f1 == 0)
      show_particle_visuals = false;
    else
      show_particle_visuals = true;
  }
  if(read1Float(line, "show_collectionBox_visuals", &f1)){
    if(f1 == 0)
      show_collectionBox_visuals = false;
    else
      show_collectionBox_visuals = true;
  }
  if(read1Float(line, "advectChoice", &f1)){
    advectChoice = (int)f1;
  }
  if(read1Float(line, "num_of_sources", &f1)){
    numOfPE = (int)f1;
    xpos = new float[numOfPE];
    ypos = new float[numOfPE];
    zpos = new float[numOfPE];
    radius = new float[numOfPE];
    rate = new float[numOfPE];
  }
  if(readSourceInfo(line, "source_info", s)){
    xpos[num] = s[0];
    ypos[num] = s[1];
    zpos[num] = s[2];
    radius[num] = s[3];
    rate[num] = s[4];
    num++;
  }
  if(read1Float(line, "release_type",&f1)){
    releaseType = (int)f1;
  }
  if(read1Float(line, "emit_method", &f1)){
    emit_method = (int)f1;
  }
  if(read1Float(line, "numBuild", &f1)){
    numBuild = (int)f1;
    xfo = new double[numBuild];
    yfo = new double[numBuild];
    zfo = new double[numBuild];
    ht = new double[numBuild];
    wti = new double[numBuild];
    lti = new double[numBuild];
  }
  if(read6Float(line, "build_param", b)){
    xfo[numb] = b[0];
    yfo[numb] = b[1];
    zfo[numb] = b[2];
    ht[numb]  = b[3];
    wti[numb] = b[4];
    lti[numb] = b[5];
    numb++;
  }
  if(read1Float(line, "pauseMode", &f1)){
    if(f1 == 0)
      pauseMode = true;
    else
      pauseMode = false;
  }
  if(read1Float(line, "calculateMeanVel", &f1)){
    if(f1 == 0)
      calculateMeanVel = false;
    else
      calculateMeanVel = true;
  }
  if(read3Float(line, "back_color", c)){
    bcolor[0] = c[0];
    bcolor[1] = c[1];
    bcolor[2] = c[2];
  }
}

bool Util::readSourceInfo(char *line, std::string settingName, float *f)
{
	std::istringstream ist(line);

	std::string w;

	ist >> w;
	if(w == settingName){	 
	  ist >> f[0];
	  ist >> f[1];
	  ist >> f[2];
	  ist >> f[3];
	  ist >> f[4];
	
	  return true;
	}
	
    return false;
}
bool Util::read3Float(char *line, std::string settingName, float *f)
{
	std::istringstream ist(line);

	std::string w;
	ist >> w;
	if(w == settingName){
		ist >> f[0];
		ist >> f[1];
		ist >> f[2];
		return true;
	}
	
    return false;
}

bool Util::read6Float(char *line, std::string settingName, float *f)
{
	std::istringstream ist(line);

	std::string w;
	ist >> w;
	if(w == settingName){
		ist >> f[0];
		ist >> f[1];
		ist >> f[2];
		ist >> f[3];
		ist >> f[4];
		ist >> f[5];
		return true;
	}
	
    return false;
}

bool Util::read1Float(char *line, std::string settingName, float *f)
{
	std::istringstream ist(line);

	std::string w;
	ist >> w;
	if(w == settingName){
		ist >> *f;
		return true;
	}
	
    return false;
}
bool Util::readComment(const char *line)
{

    if(strlen(line)==0)
    {
        return true;
    }

    int i=0;
    while(line[i] == ' ')
        i++;

    if(line[i] == '#')
    {
        return true;
    }

    return false;

}
bool Util::read1String(const char *line, char *settingName, std::string *s)
{
	std::istringstream ist(line);

	std::string w;
	ist >> w;
	if(w == settingName){
		ist >> *s;
		return true;
	}

    return false;
}
