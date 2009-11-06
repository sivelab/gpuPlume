#ifndef __QUICDATAFILE_QPSOURCE_H__
#define __QUICDATAFILE_QPSOURCE_H__ 1

#include <fstream>
#include "QUICDataFile.h"

class qpSource : public quicDataFile
{
public:
  qpSource() {}
  ~qpSource() {}

  bool readQUICFile(const std::string &filename);
  bool writeQUICFile(const std::string &filename);

!QUIC 5.51
2			!Number of sources
4			!Number of source nodes
!Start of source number 1
Source 1 !source name
1			!Source strength units (1 = g, 2 = g/s, 3 = L,4 = L/s)
15			!Source Strength
2			!Source Density (kg/m^3) [Only used for Volume based source strengths]
2			!Release Type: =1 for instantaneous;=2 for continuous; =3 for discrete continous
0			!Source start time (s)
0			!Source duration (s)
2			!Source geometry (1 = spherical shell, 2 = line, 3 = cylinder, 4 = Explosive,5 = Area, 6 = Moving Point, 7 = spherical volume, 8 = Submunitions)
2			!Number of data points
!x (m)   y (m)   z (m)
30.000000 45.000000 0.500000
60.000000 45.000000 0.500000
!End of source number 1
!Start of source number 2
Source 2 !source name
1			!Source strength units (1 = g, 2 = g/s, 3 = L,4 = L/s)
15			!Source Strength
2			!Source Density (kg/m^3) [Only used for Volume based source strengths]
2			!Release Type: =1 for instantaneous;=2 for continuous; =3 for discrete continous
0			!Source start time (s)
0			!Source duration (s)
2			!Source geometry (1 = spherical shell, 2 = line, 3 = cylinder, 4 = Explosive,5 = Area, 6 = Moving Point, 7 = spherical volume, 8 = Submunitions)
2			!Number of data points
!x (m)   y (m)   z (m)
45.000000 60.000000 0.500000
45.000000 30.000000 0.500000
!End of source number 2


private:
};

#endif // #define __QUICDATA_QPSOURCE_H__ 1
