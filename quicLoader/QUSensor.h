#ifndef __QUICDATAFILE_QUSENSOR_H__
#define __QUICDATAFILE_QUSENSOR_H__ 1

#include <cstdlib>
#include <fstream>
#include "QUICDataFile.h"
#include "legacyFileParser.h"

// //////////////////////////////////////////////////////////////////
// 
// Class for holding the QU_metparams.inp file
// 
// //////////////////////////////////////////////////////////////////
class quSensorParams : public quicDataFile
{
public:

  enum ProfileType {
    LOG = 1,
    EXP = 2,
    CANOPY = 3,
    DATAPT = 4   
  };

  std::string siteName;

  int xCoord;
  int yCoord;
  
  float decimalTime;
  ProfileType boundaryLayerFlag;

  float exponential;
  float Zo;
  float recipMoninObukhovLen;
  float canopyHeight;
  float attenuationCoef;

  float siteExponential;

  float height, speed, direction;
  
  quSensorParams() {}
  ~quSensorParams() {}

  bool readQUICFile(const std::string &filename);
  bool writeQUICFile(const std::string &filename);

private:
};

#endif // #define __QUICDATA_QUSENSOR_H__ 1
