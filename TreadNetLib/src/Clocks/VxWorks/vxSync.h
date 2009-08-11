#ifndef VXSYNC_H
#define VXSYNC_H

#include <vxWorks.h>
#include <sysLib.h>
#include <stdio.h>

class vxSync 
{
 public:
  vxSync( unsigned int hz = 32);
  ~vxSync();
  
  void block();

 protected:
  static int AnnounceCallback(int clientData);
  
 private:
  unsigned int m_hz;
  unsigned int m_ticktock;
};

inline void vxSync::block()
{
  UINT temp = m_ticktock;
  while (temp == m_ticktock) {}
}

#endif // VXSYNC_H
