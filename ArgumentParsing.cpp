/*
 *  ArgumentParsing.cpp
 *  NETCODE
 *
 *  Created by Pete Willemsen on 10/6/09.
 *  Copyright 2009 Department of Computer Science, University of Minnesota-Duluth. All rights reserved.
 *
 */

#include "ArgumentParsing.h"

ArgumentParsing::ArgumentParsing()
{
}

ArgumentParsing::ArgumentParsing(int argc, char *argv[])
{
  process(argc, argv);
}

ArgumentParsing::~ArgumentParsing()
{
}

void ArgumentParsing::reg(const std::string& argName, char shortArgName, int has_argument, bool required)
{
  ModifiedOption nextArg;

  // the option structure uses C-style character strings so get our "string" into that form.
  nextArg.optParams.name = (const char *)malloc(argName.length() + 1);
  strcpy((char *)nextArg.optParams.name, argName.c_str());

  nextArg.optParams.has_arg = has_argument;
  nextArg.optParams.flag = 0;
  nextArg.optParams.val = shortArgName;

  nextArg.isSet = false;
  nextArg.optionalArgument = "";

  m_ArgVector.push_back(nextArg);
}

bool ArgumentParsing::isSet(const std::string& argName)
{
  for (unsigned int i=0; i<m_ArgVector.size(); ++i)
  {
    if (argName.compare(m_ArgVector[i].optParams.name) == 0)
      return m_ArgVector[i].isSet;
  }
  return false;
}

bool ArgumentParsing::isSet(const std::string& argName, std::string &argValue)
{
  for (unsigned int i=0; i<m_ArgVector.size(); ++i)
    if ((argName.compare(m_ArgVector[i].optParams.name) == 0) && (m_ArgVector[i].isSet))
      {
	argValue = m_ArgVector[i].optionalArgument;
	return true;
      }

  argValue = "";
  return false;
}

int ArgumentParsing::process(int argc, char *argv[])
{
  // convert the vector into a temporary option struct needed for getopt_long
  option *getoptOptions = new option[m_ArgVector.size()];

  for (unsigned int i=0; i<m_ArgVector.size(); ++i)
    {
      getoptOptions[i].name = (char *)malloc(strlen(m_ArgVector[i].optParams.name) + 1);
      strcpy((char *)getoptOptions[i].name, m_ArgVector[i].optParams.name);
      getoptOptions[i].has_arg = m_ArgVector[i].optParams.has_arg;
      getoptOptions[i].flag = m_ArgVector[i].optParams.flag;
      getoptOptions[i].val = m_ArgVector[i].optParams.val;
    }

  int c;
  while (1)
    {
      // which option index are we on?
      int this_option_optind = optind ? optind : 1;
      int option_index = 0;

      c = getopt_long(argc, argv, "", getoptOptions, &option_index);

      // when the getopt functions return a -1, there are no more
      // arguments to process.
      if (c == -1)
	break;

      int argIdx = 0;
      bool found = false;
      while (!found && argIdx < m_ArgVector.size() && c != '?')
	{
	  if (c == (int)m_ArgVector[argIdx].optParams.val)
	    {
	      // std::cout << "found option: " << m_ArgVector[argIdx].optParams.name << std::endl;
	      m_ArgVector[argIdx].isSet = true;
	      if (m_ArgVector[argIdx].optParams.has_arg != no_argument)
		{
		  m_ArgVector[argIdx].optionalArgument = optarg;	      
		  // std::cout << "\t optarg = " << optarg << std::endl;
		}
	      
	      found = true;
	    }

	  argIdx++;
	}

      if (!found)
	{
	  std::cerr << "?? getopt returned character code: " << c << std::endl;
	}
    }

#if 0
  // look over the non-option elements
  if (optind < argc) {
    printf("non-option ARGV-elements: ");
    while (optind < argc)
      printf("%s ", argv[optind++]);
    printf("\n");
  }
#endif

  // deallocate the memory we created to make this happen
  for (unsigned int i=0; i<m_ArgVector.size(); ++i)
  {
      free((char *)(getoptOptions[i].name));
  }
  
  delete [] getoptOptions;

  return 1;
}
