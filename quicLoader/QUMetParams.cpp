#include "QUMetParams.h"

bool quMetParams::readQUICFile(const std::string &filename)
{
  std::cout << "\tParsing QU_metparams.inp file: " << filename << std::endl;

  std::ifstream quicFile(filename.c_str(), std::ifstream::in);
  if(!quicFile.is_open())
    {
      std::cerr << "quicLoader could not open :: " << filename << "." << std::endl;
      exit(EXIT_FAILURE);
    }
		
  std::string line;
  std::istringstream ss(line, std::stringstream::in);

  // first thing in these files is now a comment about the version...
  getline(quicFile, line);

  std::string quicComment, quicVersion;
  ss.str(line);
  ss >> quicComment >> quicVersion;
  // std::cout << "quicComment: " << quicComment << ", quicVersion: " << quicVersion << std::endl;
  if (atof(quicVersion.c_str()) != 5.72)
    {
      std::cerr << "Error!  Exiting!  Only QUIC version 5.72 supported!  Read file with Version " << quicVersion << std::endl;
      exit(EXIT_FAILURE);
    }
  ss.clear();

  // !Met input flag (0=QUIC,1=ITT MM5,2=HOTMAC)
  getline(quicFile, line);
  ss.str(line);

  // This little piece of code is here because strange
  // end-of-line/carriage return characters are showing up after the
  // version number.  In particular, 0x0a30.  For now, I've dealt with
  // this by just clearing the stringstream's buffer.  Seems to remove
  // issues with accumulating strings.
  int mit = -1;
  ss >> mit;
  ss.clear();

  if (mit == quMetParams::QUIC)
    metInputFlag = quMetParams::QUIC;
  else if (mit == quMetParams::ITT_MM5)
    metInputFlag = quMetParams::ITT_MM5;
  else if (mit == quMetParams::HOTMAC)
    metInputFlag = quMetParams::HOTMAC;
  else 
    {
      std::cout << "quicLoader: unknown Met Input Flag type provided: " << mit << std::endl;
      exit(EXIT_FAILURE);
    }
		
  // !Number of measuring sites
  getline(quicFile, line);
  ss.str(line);
  ss >> numMeasuringSites;
  ss.clear();


  // !Maximum size of data points profiles
  getline(quicFile, line);
  ss.str(line);
  ss >> maxSizeProfiles;
  ss.clear();


		
  // !Site Name 
  getline(quicFile, line);
  ss.str(line);
  ss >> siteName;
		
  // Need to skip over this non-standard format... which has the !File name comment prior to 
  // the actual filename... argh...
  //
  // !File name
  getline(quicFile, line);
		
  // the actual !file name 
  getline(quicFile, line);
  ss.str(line);
  ss >> fileName;

  quicFile.close();
  return true;
}

bool quMetParams::writeQUICFile(const std::string &filename)
{
  std::ofstream qufile;
  qufile.open(filename.c_str());
  qufile << "!QUIC 5.72" << std::endl;

  if (qufile.is_open())
    {
      qufile << metInputFlag << "\t\t\t!Met input flag (0=QUIC,1=ITT MM5,2=HOTMAC)" << std::endl;
      qufile << numMeasuringSites << "\t\t\t!Number of measuring sites" << std::endl;
      qufile << maxSizeProfiles << "\t\t\t!Maximum size of data points profiles" << std::endl;
      qufile << siteName << "\t\t\t!Site Name " << std::endl;
      qufile << "!File name" << std::endl;
      qufile << fileName << std::endl;

      return true;
    }

  return false;
}
