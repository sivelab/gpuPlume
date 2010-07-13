#include "QUSensor.h"

bool quSensorParams::readQUICFile(const std::string &filename)
{
  std::cout << "\tParsing sensor file: " << filename << std::endl;

  std::ifstream quicFile(filename.c_str(), std::ifstream::in);
  if(!quicFile.is_open())
    {
      std::cerr << "quicLoader could not open :: " << filename << "." << std::endl;
      exit(EXIT_FAILURE);
    }
		
  std::string line;
  std::stringstream ss(line, std::stringstream::in | std::stringstream::out);

  getline(quicFile, line);
  ss.str(line);
  ss >> siteName;
		
  getline(quicFile, line);
  ss.str(line);
  ss >> xCoord;

  getline(quicFile, line);
  ss.str(line);
  ss >> yCoord;

  getline(quicFile, line);
  ss.str(line);
  ss >> decimalTime;

  getline(quicFile, line);
  ss.str(line);
  int profileType;
  ss >> profileType;

  boundaryLayerFlag = (ProfileType)profileType;
  if (boundaryLayerFlag == LOG)
    {
      getline(quicFile, line);
      ss.str(line);
      ss >> Zo;
      
      getline(quicFile, line);
      ss.str(line);
      ss >> recipMoninObukhovLen;
    }
  else if (boundaryLayerFlag == EXP)
    {
      getline(quicFile, line);
      ss.str(line);
      ss >> exponential;
    }
  else if (boundaryLayerFlag == CANOPY)
    {
      getline(quicFile, line);
      ss.str(line);
      ss >> Zo;
      
      getline(quicFile, line);
      ss.str(line);
      ss >> recipMoninObukhovLen;

      getline(quicFile, line);
      ss.str(line);
      ss >> canopyHeight;
      
      getline(quicFile, line);
      ss.str(line);
      ss >> attenuationCoef;
    }
  else if (boundaryLayerFlag == DATAPT)
    {
      std::cerr << "Do not yet support Discrete Data Points Wind Profiles yet!!!!" << std::endl;
      exit(EXIT_FAILURE);
    }
  else
    {
      std::cerr << "Unknown Wind Profile type: " << boundaryLayerFlag << "!!! Exiting!" << std::endl;
      exit(EXIT_FAILURE);
    }

  // step over the height, speed, direction header
  getline(quicFile, line);
  ss.str(line);

  getline(quicFile, line);
  ss.str(line);
  ss >> height >> speed >> direction;

  quicFile.close();
  return true;
}

bool quSensorParams::writeQUICFile(const std::string &filename)
{
  std::ofstream qufile;
  qufile.open(filename.c_str());

  if (qufile.is_open())
    {
      qufile << siteName << "\t\t\t!Site Name " << std::endl;
      qufile << xCoord << "\t\t\t!X coordinate (meters)" << std::endl;
      qufile << yCoord << "\t\t\t!Y coordinate (meters)" << std::endl;
      qufile << decimalTime << "\t\t\t!Decimal time (military time i.e. 0130 = 1.5)" << std::endl;
      qufile << boundaryLayerFlag << "\t\t\t!site boundary layer flag (1 = log, 2 = exp, 3 = urban canopy, 4 = discrete data points)" << std::endl;

      if (boundaryLayerFlag == LOG)
	{
	  qufile << Zo << "\t\t\t!site zo" << std::endl;
	  qufile << recipMoninObukhovLen << "\t\t\t!reciprocal Monin-Obukhov Length (1/m)" << std::endl;
	}
      else if (boundaryLayerFlag == EXP)
	{
	  qufile << exponential << "\t\t\t!site exponential" << std::endl;
	}
      else if (boundaryLayerFlag == CANOPY)
	{
	  qufile << Zo << "\t\t\t!site zo" << std::endl;
	  qufile << recipMoninObukhovLen << "\t\t\t!reciprocal Monin-Obukhov Length (1/m)" << std::endl;
	  qufile << canopyHeight << "\t\t\t!Canopy height (m)" << std::endl;
	  qufile << attenuationCoef << "\t\t\t!attenuation coefficient" << std::endl;
	}
      else if (boundaryLayerFlag == DATAPT)
	{
	  std::cerr << "Do not yet support Discrete Data Points Wind Profiles yet!!!!" << std::endl;
	  exit(EXIT_FAILURE);
	}
      else
	{
	  std::cerr << "Unknown Wind Profile type: " << boundaryLayerFlag << "!!! Exiting!" << std::endl;
	  exit(EXIT_FAILURE);
	}

      qufile << "!Height (m),Speed	(m/s), Direction (deg relative to true N)" << std::endl;
      qufile << height << " " << speed << " " << direction << std::endl;

      return true;
    }

  return false;
}
