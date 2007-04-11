#ifndef __UTIL_H__
#define __UTIL_H__

#include <iostream>
#include <fstream>

class PlumeControl;


class Util{

 public:
  
  Util(PlumeControl*);
  void readInput(std::string);

 private:

  void parseLine(char*);
  bool read1Float(char*,std::string,float*);
  bool read6Float(char*,std::string,float*);
  bool readComment(const char*);
  bool read1String(const char*,char*,std::string*);
  
  PlumeControl* plume;

};

#endif //__UTIL_H__
