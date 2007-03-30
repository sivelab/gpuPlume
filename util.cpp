#include "util.h"
#include <stdio.h>
#include "plumeControl.h"

Util::Util(PlumeControl* p){
  plume = p;
}

void Util::readInput(char* file){
  FILE *in = fopen(file,"r");

  
  char line[1024];
  while( fgets(line, 1024, in) != NULL )
  {
      if( line[ strlen(line)-1 ] == '\n' )
          line[ strlen(line)-1 ] = '\0';

       parseLine(line);
  }
    
  fclose(in);

}
void Util::parseLine(char* line){
  char s1[1024];
  float f1;

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
  if(read1String(line, "output_file", s1)){
    plume->output_file = std::string(s1);
  }


}
bool Util::read1Float(const char *line, char *settingName, float *f)
{
    char readSettingName[256];
    char readSettingVal[512];
    readSettingName[0] = readSettingVal[0] = 0;

    if(sscanf(line, "%s %[^\n]\n", readSettingName, readSettingVal) != 2)
        return false;

    if(strcmp(readSettingName, settingName) == 0)
    {
        if(sscanf(readSettingVal, "%f", f) != 1)
        {
           exit(1);
        }
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
bool Util::read1String(const char *line, char *settingName, char *s)
{
    char readSettingName[256];
    char readSettingVal[512];
    readSettingName[0] = readSettingVal[0] = 0;

    if(sscanf(line, "%s %[^\n]\n", readSettingName, readSettingVal) != 2)
        return false;

    if(strcmp(readSettingName, settingName) == 0)
    {
        if(sscanf(readSettingVal, "%s", s) != 1)
        {
            exit(1);
        }
        return true;
    }
    return false;
}
