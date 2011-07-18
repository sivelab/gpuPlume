/*
 *  ArgumentParsing.h
 *  NETCODE
 *
 *  Created by Pete Willemsen on 10/6/09.
 *  Copyright 2009 Department of Computer Science, University of Minnesota-Duluth. All rights reserved.
 *
 */

#ifndef __ARGUMENT_PARSING_H__
#define __ARGUMENT_PARSING_H__ 1

#include <iostream>
#include <cstdlib>
#include <cstring>
#include <vector>

#ifdef WIN32
#include "getopt_win32.h"
#else
#include <getopt.h>
#endif

class ArgumentParsing
{
public:
  ArgumentParsing();
  ArgumentParsing(int argc, char *argv[]);
  ~ArgumentParsing();

  void reg(const std::string& argName, char shortArgName, int has_argument, bool required=false);
  int processCommandLineArgs(int argc, char *argv[]) { return process(argc, argv); }

  bool isSet(const std::string& argName);
  bool isSet(const std::string& argName, std::string &argValue);

protected:
  int process(int argc, char *argv[]);
	
private:

  // getopt_long structure variables
  // struct option {
  //    const char *name;
  //    int         has_arg;
  //    int        *flag;
  //    int         val;
  // };

  struct ModifiedOption 
  {
	bool isSet;
    option optParams;
    std::string optionalArgument;
  };
   
  std::vector<ModifiedOption> m_ArgVector;
};

#endif // __ARGUMENT_PARSING_H__ 1
