#ifndef __QUICDATAFILE_QUBUILDINGS_H__
#define __QUICDATAFILE_QUBUILDINGS_H__ 1

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <vector>

#include "QUICDataFile.h"

// //////////////////////////////////////////////////////////////////
// 
// Class for holding the QU_buildings.inp file
// 
// //////////////////////////////////////////////////////////////////
class quBuildings : public quicDataFile
{
public:
  quBuildings() {}
  ~quBuildings() {}

  bool readQUICFile(const std::string &filename);
  bool writeQUICFile(const std::string &filename);

  int x_subdomain_sw;
  int y_subdomain_sw;
  int x_subdomain_ne;
  int y_subdomain_ne;

  float zo;

  struct buildingData
  {
    int bldNum;
    int group;
    int type;
    float height;
    float width;
    float length;
    float xfo;
    float yfo;
    float zfo;
    float gamma;
    float supplementalData;
  };

  std::vector<buildingData> buildings;
  
private:
};

#endif // #define __QUICDATAFILE_QUBUILDINGS_H__
