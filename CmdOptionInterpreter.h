#ifndef CMD_OPTION_INTERPRETER_H
#define CMD_OPTION_INTERPRETER_H

/**
 *  CmdOptionInterpreter is responcible for parsing and setting 
 *  options that are passed in via the command line when the 
 *  program starts.
 */

#include <cassert>
#include <string>
#include <cstdlib>
#include <cmath>

#include "util/ArgumentParsing.h"
#include "util.h"

class CmdOptionInterpreter {

 public:
  
  /**
   * Default Constructor
   */
  CmdOptionInterpreter();
  
  /**
   *  Constructor that takes in a reference 
   *  to a Util object.
   */
  CmdOptionInterpreter(sivelab::ArgumentParsing *argParser, Util * util);
  
  /**
   *  Destructor.
   */
  ~CmdOptionInterpreter();

  /**
   * parse takes in the input from the console when the program
   * was launched and parses each option. It then sets those
   * passed options within the given Util object.
   *
   * This method should be run just after the Input.txt file
   * has been processed. This way these settings will override
   * any Input.txt file.
   */
  void parse();
  
 private:

  Util *m_utilPtr;
  sivelab::ArgumentParsing *m_argParser;

};

#endif // CMD_OPTION_INTERPRETER_H
