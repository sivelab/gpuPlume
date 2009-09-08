
/**
 *  CmdOptionInterpreter is responcible for parsing and setting 
 *  options that are passed in via the command line when the 
 *  program starts.
 */

#include "CmdOptionInterpreter.h"

CmdOptionInterpreter::CmdOptionInterpreter() {
  data = NULL;
}

CmdOptionInterpreter::CmdOptionInterpreter(Util * util) {
  data = util;
}

CmdOptionInterpreter::~CmdOptionInterpreter() {
  
}

void CmdOptionInterpreter::parse(int argc, char ** argv) {
  
  // Check to make sure that the correct constructor was
  // used and that data is set.
  if(data == NULL) {
    return;
  }
  
  // Loop through the remaining options that were passed.
  // Note that we assume that the first option is a .proj
  // file and that it has already been read.
  for(int i = 2; i < argc; i++) {
    std::string option = std::string(argv[i]);
    
    if(option.compare("--fullscreen") == 0) {
      data->fullscreen = true;
    } else if (option.compare("--networkmode") == 0) {
      std::string variable = std::string(argv[i+1]);
      data->network_mode = atoi(variable.c_str());
    } else if (option.compare("--viewingmode") == 0) {
      std::string variable = std::string(argv[i+1]);
      data->viewing_mode = atoi(variable.c_str());
      i++;
    } else if (option.compare("--treadportview") == 0) {
      std::string variable = std::string(argv[i+1]);
      data->treadport_view = variable.c_str()[0];
    } else if (option.compare("--dynamicTreadportFrustum") == 0) {
      data->static_treadport_frustum = 0;
    } else if (option.compare("--sunAzimuth") == 0) {
      std::string variable = std::string(argv[i+1]);
      data->sun_azimuth = atoi(variable.c_str());
    } else if (option.compare("--sunAltitude") == 0) {
      std::string variable = std::string(argv[i+1]);
      data->sun_altitude = atoi(variable.c_str());
    } else if (option.compare("--onlyCalcShadows") == 0) {
      data->onlyCalcShadows = true;
    }
  }
  
  return;
}
