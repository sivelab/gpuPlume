#include "vxSync.h"

vxSync::vxSync( unsigned int hz)
  : m_hz( hz )
{
  m_ticktock = 0;
  char emsg[128];

  // set rate of aux clk
  if (sysAuxClkRateSet(m_hz) == ERROR) {
    m_hz = sysAuxClkRateGet();
    sprintf(emsg, "*** sysAuxClkRateSet failed - currently %i\n", m_hz);
    perror(emsg);
  }
  
  // connect aux clk to function
  if (sysAuxClkConnect((FUNCPTR) vxSync::AnnounceCallback,
		       (int) this) == ERROR) {
    sprintf(emsg, "**** Error: sysAuxClkConnect failed.\n");
    perror(emsg);
  }

  // turn clock interrupts on
  sysAuxClkEnable();
}

vxSync::~vxSync()
{
  sysAuxClkDisable();
}

int vxSync::AnnounceCallback(int clientData) {
  vxSync *my_sync = (vxSync *)clientData;

  my_sync->m_ticktock = !my_sync->m_ticktock;
}
