#ifndef WINSYNC_H
#define WINSYNC_H

#include <windows.h>
#include <mmsystem.h>
#include <stdio.h>

class winSync 
{
 public:
  winSync( UINT hz = 32, UINT res = 1 );
  ~winSync();
  
  void block();
  
 protected:
  static void CALLBACK clockTimerCB(UINT uID, UINT uMsg, DWORD dwUser,
				    DWORD dw1, DWORD dw2);
  
 private:
  UINT m_timerID;
  UINT m_hz;
  UINT m_res;
  double m_true_hz;
  UINT m_ticktock;
};

inline void winSync::block()
{
  UINT temp = m_ticktock;
  while (temp == m_ticktock) {}
}

#endif // WIN_SYNC_H
