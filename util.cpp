#include "util.h"
#include <stdio.h>
#include "gpuPlume.h"

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
    sigU = 2.0f*f1;
    sigV = 2.0f*f1;
    sigW = 1.3f*f1;
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
    xfo = new float[numBuild];
    yfo = new float[numBuild];
    zfo = new float[numBuild];
    ht = new float[numBuild];
    wti = new float[numBuild];
    lti = new float[numBuild];
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
