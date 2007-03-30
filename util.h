#ifndef __UTIL_H__
#define __UTIL_H__

#include <iostream>
#include <fstream>
// #include "plumeControl.h"

class PlumeControl;

class Util{

 public:
  
  Util(PlumeControl*);
  void readInput(char*);

 private:

  void parseLine(char*);
  bool read1Float(const char*,char*,float*);
  bool readComment(const char*);
  bool read1String(const char*,char*,char*);
  
  PlumeControl* plume;

};

#endif //__UTIL_H__
