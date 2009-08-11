#include "winSync.h"

winSync::winSync( UINT hz, UINT res )
  : m_hz( hz ), m_res( res )
{
  char emsg[128];
  
  if (m_res < 0 || m_res > 1000000) {
    sprintf(emsg, "*** Error : specified resolution is out of range\n");
    perror(emsg);
  }

  // convert rate to millesecs
  UINT sleepTime = (UINT) (1000.0 / (double) m_hz);
  m_true_hz = (double)1000.0/sleepTime;

  m_timerID = timeSetEvent(sleepTime, m_res,
			   clockTimerCB, (ULONG)this,
			   TIME_PERIODIC | TIME_CALLBACK_FUNCTION);

  if (m_timerID == NULL) {
    sprintf(emsg, "*** Error : Can't initialize multimedia timer\n");
    perror(emsg);
  }
}

winSync::~winSync()
{
  timeKillEvent(m_timerID);
}

void CALLBACK winSync::clockTimerCB(UINT timerID, UINT uMsg,
				    DWORD self, DWORD dw1, DWORD dw2) {
  winSync *my_sync = (winSync *)self;

  my_sync->m_ticktock = !my_sync->m_ticktock;
}
