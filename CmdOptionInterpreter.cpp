
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

CmdOptionInterpreter::CmdOptionInterpreter(sivelab::ArgumentParsing *argParser, Util *util) 
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

  if (m_argParser->isSet("version"))
    {
      std::cout << "gpuPlume" << std::endl;
      // std::cout << "Version: " << GPUPLUME_VERSION << std::endl;
    }

  if (m_argParser->isSet("fullscreen"))
    m_utilPtr->fullscreen = true;

  int tmpInt;
  if (m_argParser->isSet("networkmode", tmpInt))
    m_utilPtr->network_mode = tmpInt;

  tmpInt = 0;
  if (m_argParser->isSet("viewingmode", tmpInt)) 
    m_utilPtr->viewing_mode = tmpInt;

  char tmpChar;
  if (m_argParser->isSet("treadportview", tmpChar))
    m_utilPtr->treadport_view = tmpChar;

  if (m_argParser->isSet("dynamicTreadportFrustum"))
    m_utilPtr->static_treadport_frustum = 0;

  float tmpFloat = 0.0;
  if (m_argParser->isSet("sunAzimuth", tmpFloat))
    m_utilPtr->sun_azimuth = tmpFloat;

  tmpFloat = 0.0;
  if (m_argParser->isSet("sunAltitude", tmpFloat))
    m_utilPtr->sun_altitude = tmpFloat;

  if (m_argParser->isSet("onlyCalcShadows"))
    m_utilPtr->onlyCalcShadows = true;

  tmpInt = 0;
  if (m_argParser->isSet("numParticles", tmpInt))
    {
      m_utilPtr->qpParamData.numParticles = tmpInt;
      m_utilPtr->twidth = (int)sqrt( static_cast<float>(m_utilPtr->qpParamData.numParticles) );
      m_utilPtr->theight = (int)sqrt( static_cast<float>(m_utilPtr->qpParamData.numParticles) );
      std::cout << "COMMAND LINE OVERRIDE: using num particles=" << m_utilPtr->qpParamData.numParticles << ". Actually using " << m_utilPtr->twidth * m_utilPtr->theight << " particles!" << std::endl;
    }

  std::string tmpString = "";
  if (m_argParser->isSet("concFile", tmpString))
    {
      m_utilPtr->output_file = tmpString;
      std::cout << "COMMAND LINE OVERRIDE: using concentration output file: " << m_utilPtr->output_file << std::endl;
    }

  tmpString = "";
  if (m_argParser->isSet("concId", tmpString))
    {
      m_utilPtr->output_id = tmpString;
      std::cout << "COMMAND LINE OVERRIDE: using concentration id = " << m_utilPtr->output_id << std::endl;
    }

  if (m_argParser->isSet("ignoreSignal"))
    m_utilPtr->ignoreSignal = true;

  if (m_argParser->isSet("offscreenRender"))
    m_utilPtr->offscreenRender = true;

  tmpInt = 0;
  if (m_argParser->isSet("problemID", tmpInt))
    m_utilPtr->problemID = tmpInt;

  tmpInt = 0;
  if (m_argParser->isSet("probInstID", tmpInt))
    m_utilPtr->problemInstanceID = tmpInt;
}
