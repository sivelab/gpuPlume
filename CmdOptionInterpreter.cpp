
/**
 *  CmdOptionInterpreter is responcible for parsing and setting 
 *  options that are passed in via the command line when the 
 *  program starts.
 */

#include "CmdOptionInterpreter.h"

CmdOptionInterpreter::CmdOptionInterpreter() 
{
  m_utilPtr = NULL;
  m_argParser = NULL;
}

CmdOptionInterpreter::CmdOptionInterpreter(ArgumentParsing *argParser, Util *util) 
{
  m_utilPtr = util;
  m_argParser = argParser;
}

CmdOptionInterpreter::~CmdOptionInterpreter() 
{
}

void CmdOptionInterpreter::parse() 
{
  // Check to make sure that the correct constructor was
  // used and that m_utilPtr is set.
  assert(m_utilPtr && m_argParser);
  
  std::string argVal = "";

  if (m_argParser->isSet("fullscreen"))
    m_utilPtr->fullscreen = true;

  if (m_argParser->isSet("networkmode", argVal))
    m_utilPtr->network_mode = atoi(argVal.c_str());

  if (m_argParser->isSet("viewingmode", argVal)) 
    m_utilPtr->viewing_mode = atoi(argVal.c_str());

  if (m_argParser->isSet("treadportview", argVal))
    m_utilPtr->treadport_view = argVal.c_str()[0];

  if (m_argParser->isSet("dynamicTreadportFrustum"))
    m_utilPtr->static_treadport_frustum = 0;

  if (m_argParser->isSet("sunAzimuth", argVal))
    m_utilPtr->sun_azimuth = atoi(argVal.c_str());

  if (m_argParser->isSet("sunAltitude", argVal))
    m_utilPtr->sun_altitude = atoi(argVal.c_str());

  if (m_argParser->isSet("onlyCalcShadows"))
    m_utilPtr->onlyCalcShadows = true;

  if (m_argParser->isSet("numParticles", argVal))
    {
      m_utilPtr->qpParamData.numParticles = atoi(argVal.c_str());
      m_utilPtr->twidth = (int)sqrt(m_utilPtr->qpParamData.numParticles);
      m_utilPtr->theight = (int)sqrt(m_utilPtr->qpParamData.numParticles);
      std::cout << "COMMAND LINE OVERRIDE: using num particles=" << m_utilPtr->qpParamData.numParticles << ". Actually using " << m_utilPtr->twidth * m_utilPtr->theight << " particles!" << std::endl;
    }

  if (m_argParser->isSet("concFile", argVal))
    {
      m_utilPtr->output_file = argVal;
      std::cout << "COMMAND LINE OVERRIDE: using concentration output file: " << m_utilPtr->output_file << std::endl;
    }

  if (m_argParser->isSet("concId", argVal))
    {
      m_utilPtr->output_id = argVal;
      std::cout << "COMMAND LINE OVERRIDE: using concentration id = " << m_utilPtr->output_id << std::endl;
    }
}
