#include "util.h"
#include <stdio.h>
#include <sstream>
#include "plumeControl.h"

Util::Util(PlumeControl* p){
	plume = p;
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
  //char s1[1024];
  float f1;
  std::string s1;

  if(readComment(line))
    return;

  if(read1Float(line, "twidth", &f1)){
	 plume->twidth = (int)f1;
  }
  if(read1Float(line, "theight", &f1)){
	 plume->theight = (int)f1;
  }
  if(read1Float(line, "time_step", &f1)){
    plume->time_step = f1;
  }
  if(read1Float(line, "testcase", &f1)){
    plume->testcase = (int)f1;
  }
  if(read1Float(line, "useRealTime", &f1)){
    if(f1 == 0)
      plume->useRealTime = false;
    else
      plume->useRealTime = true;
  }
  if(read1String(line, "output_file", &s1)){
	  plume->output_file = s1;
  }
  if(read1Float(line, "duration", &f1)){
    plume->duration = (double)f1;
  }

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
