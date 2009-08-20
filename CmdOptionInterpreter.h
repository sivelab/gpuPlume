#ifndef CMD_OPTION_INTERPRETER_H
#define CMD_OPTION_INTERPRETER_H

/**
 *  CmdOptionInterpreter is responcible for parsing and setting 
 *  options that are passed in via the command line when the 
 *  program starts.
 */

#include "util.h"

#include <string>

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
  CmdOptionInterpreter(Util * util);
  
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
  void parse(int argc, char ** argv);
  
 private:

  Util * data;

};

#endif // CMD_OPTION_INTERPRETER_H
